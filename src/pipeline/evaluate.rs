/// Main assembly evaluation pipeline
/// Orchestrates: FASTA validation → read mapping → SV detection → error merging → QV computation

use anyhow::Result;
use log::{info, warn};
use std::process::Command;
use rayon::prelude::*;
use crate::{static_analysis, detect, merge, base_error, plot, utils};

#[derive(Debug, Clone)]
pub struct EvaluateConfig {
    pub contig: Vec<String>,
    pub read: Vec<String>,
    pub datatype: String,
    pub outpath: String,
    pub reference: Option<String>,
    pub thread: usize,
    pub read_coverage: usize,
    pub min_depth: Option<usize>,
    pub min_contig_length: usize,
    pub min_contig_length_assemblyerror: usize,
    pub min_assembly_error_size: usize,
    pub max_assembly_error_size: usize,
    pub noplot: bool,
    pub skip_read_mapping: bool,
    pub skip_structural_error: bool,
    pub skip_structural_error_detect: bool,
    pub skip_base_error: bool,
    pub skip_base_error_detect: bool,
}

/// Main entry point for assembly evaluation
pub fn run(config: EvaluateConfig) -> Result<()> {
    let start_time = std::time::Instant::now();

    // Validate inputs
    info!("Validating input files");
    utils::validate_input_files(&config.contig)?;
    utils::validate_input_files(&config.read)?;
    utils::validate_datatype(&config.datatype)?;

    // Normalize output path
    let outpath = utils::normalize_path(&config.outpath);
    utils::ensure_output_dir(&outpath)?;

    info!("Output directory: {}", outpath);

    // Initialize logger
    let log_file = format!("{}Inspector.log", outpath);
    let mut log = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&log_file)?;

    use std::io::Write;
    writeln!(log, "Inspector starting... {}", chrono::Local::now().format("%d/%m/%Y %H:%M:%S"))?;
    writeln!(log, "Start Assembly evaluation with contigs: {:?}", config.contig)?;

    // Phase 2: Static analysis (FASTA validation)
    info!("Phase 2: FASTA validation and statistics");
    let fasta_stats = static_analysis::simple_fasta(
        &config.contig,
        &outpath,
        config.min_contig_length,
        config.min_contig_length_assemblyerror,
    )?;

    info!("Total contig length: {}", utils::format_size(fasta_stats.total_length));
    info!("N50: {}", fasta_stats.n50);
    info!("Large contigs: {}", fasta_stats.chromosomes_large.len());

    let t1 = std::time::Instant::now();
    writeln!(log, "TIME: Before read mapping {}", (t1 - start_time).as_secs_f64())?;

    // Phase 3a: Read mapping (via minimap2 + samtools)
    if !config.skip_read_mapping {
        info!("Phase 3a: Mapping reads to contigs");
        map_reads(&config, &outpath, fasta_stats.total_length, &mut log)?;
    }

    let t2 = std::time::Instant::now();
    writeln!(log, "TIME: Read Alignment: {}", (t2 - t1).as_secs_f64())?;

    // Phase 3b/4: SV detection and merging
    if !config.skip_structural_error_detect {
        info!("Phase 3b & 4: SV detection and merging");
        
        // Create workspace
        std::fs::create_dir_all(format!("{}map_depth/", outpath))?;
        
        if !config.skip_structural_error {
            std::fs::create_dir_all(format!("{}debreak_workspace/", outpath))?;

            // Parallel SV detection over large contigs (each writes to separate per-chrom file)
            fasta_stats.chromosomes_large.par_iter()
                .map(|chrom| detect::detect_sortbam(&outpath, config.min_assembly_error_size, config.max_assembly_error_size, chrom))
                .collect::<Result<Vec<_>>>()?;

            // Parallel depth collection for small contigs
            let small_chroms: Vec<_> = fasta_stats.chromosomes_map.keys()
                .filter(|c| !fasta_stats.chromosomes_large.contains(c))
                .collect();
            small_chroms.par_iter()
                .map(|chrom| detect::detect_sortbam_nosv(&outpath, chrom, "small"))
                .collect::<Result<Vec<_>>>()?;
        }
    }

    let coverage = static_analysis::mapping_info_contig(
        &outpath,
        &fasta_stats.chromosomes_large,
        &fasta_stats.chromosomes_map.keys()
            .filter(|c| !fasta_stats.chromosomes_large.contains(c))
            .cloned()
            .collect::<Vec<_>>(),
        fasta_stats.total_length,
        fasta_stats.total_length_large,
    )?;

    let min_support = std::cmp::max(1, coverage / 10);
    info!("Average coverage: {}", coverage);
    info!("Min support: {}", min_support);

    let t3 = std::time::Instant::now();
    writeln!(log, "TIME: Structural error signal detection: {}", (t3 - t2).as_secs_f64())?;

    // Phase 5: SV clustering and filtering
    let ae_len_structural_error;
    if !config.skip_structural_error {
        std::fs::create_dir_all(format!("{}ae_merge_workspace", outpath))?;
        
        // Parallel SV clustering — each contig reads/writes to separate files
        fasta_stats.chromosomes_large.par_iter()
            .map(|chrom| {
                let contig_length = fasta_stats.contig_lengths.get(chrom).copied().unwrap_or(0);
                merge::cluster(&outpath, chrom, contig_length, min_support, coverage * 2)?;
                merge::cluster_insertions(&outpath, chrom, contig_length, min_support, coverage * 2, "ins")?;
                merge::cluster_insertions(&outpath, chrom, contig_length, min_support, coverage * 2, "inv")?;
                Ok(())
            })
            .collect::<Result<Vec<_>>>()?;
        
        static_analysis::assembly_info_cluster(&outpath, config.min_assembly_error_size, config.max_assembly_error_size)?;
        merge::genotype(coverage, &outpath)?;
        ae_len_structural_error = merge::filter_errors(coverage, &outpath, config.min_assembly_error_size, &config.datatype)?;
    } else {
        ae_len_structural_error = 0;
    }

    let t4 = std::time::Instant::now();
    writeln!(log, "TIME: Structural error clustering : {}", (t4 - t3).as_secs_f64())?;

    // Phase 6: Base error detection
    let ae_len_base_error;
    if !config.skip_base_error {
        info!("Phase 6: Base error detection");
        
        if !config.skip_base_error_detect {
            std::fs::create_dir_all(format!("{}base_error_workspace", outpath))?;

            // Parallel pileup SNV detection per contig (each writes to separate baseerror_{chrom}.bed)
            let all_chroms: Vec<_> = fasta_stats.chromosomes_map.keys().collect();
            all_chroms.par_iter()
                .map(|chrom| base_error::get_snv(
                    &outpath,
                    chrom,
                    (coverage * 2) / 5,
                    coverage * 2,
                    config.min_depth,
                ))
                .collect::<Result<Vec<_>>>()?;
        }
        
        ae_len_base_error = base_error::count_base_errors(
            &outpath,
            fasta_stats.total_length,
            &config.datatype,
            coverage,
        )?;
    } else {
        ae_len_base_error = 0;
    }

    let t5 = std::time::Instant::now();
    writeln!(log, "TIME: Small-scale error detection: {}", (t5 - t4).as_secs_f64())?;

    // Compute long_read_QV — always, regardless of error count
    {
        use std::io::Write;
        let total_errors = ae_len_structural_error + ae_len_base_error;
        let valid_bases = fasta_stats.total_length;

        let long_read_qv = if total_errors == 0 || valid_bases == 0 {
            f64::INFINITY
        } else {
            -10.0 * (total_errors as f64 / valid_bases as f64).log10()
        };

        // Write long_read_QV to summary_statistics (already created in Phase 2; append new entries)
        let mut summary = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(format!("{}summary_statistics", outpath))?;
        writeln!(summary, "Assembly error size\t{}", total_errors)?;
        writeln!(summary, "Structural error\t{}", ae_len_structural_error)?;
        writeln!(summary, "Base error\t{}", ae_len_base_error)?;
        if long_read_qv.is_infinite() {
            writeln!(summary, "long_read_QV\tInf")?;
        } else {
            writeln!(summary, "long_read_QV\t{:.2}", long_read_qv)?;
        }

        // Log long_read_QV prominently
        if long_read_qv.is_infinite() {
            info!("Assembly error size: {} (structural: {}, base: {})", total_errors, ae_len_structural_error, ae_len_base_error);
            info!("long_read_QV: Inf (no errors detected)");
        } else {
            info!("Assembly error size: {} (structural: {}, base: {})", total_errors, ae_len_structural_error, ae_len_base_error);
            info!("long_read_QV: {:.2}", long_read_qv);
        }

        writeln!(log, "Assembly error size: {}", total_errors)?;
        if long_read_qv.is_infinite() {
            writeln!(log, "long_read_QV: Inf")?;
        } else {
            writeln!(log, "long_read_QV: {:.2}", long_read_qv)?;
        }
    }

    let t6 = std::time::Instant::now();
    writeln!(log, "TIME: long_read_QV calculation: {}", (t6 - t5).as_secs_f64())?;

    // Phase 7: Plotting
    if !config.noplot {
        info!("Phase 7: Generating plots - skipping");
        plot::plot_n50(&outpath, config.min_contig_length)?;
        plot::plot_dotplot(&outpath)?;
    }

    let total_time = std::time::Instant::now();
    writeln!(log, "Inspector finished. Total time: {}", (total_time - start_time).as_secs_f64())?;

    info!("Assembly evaluation complete in {:.2}s", (total_time - start_time).as_secs_f64());
    Ok(())
}

/// Wrapper for read-to-contig alignment using minimap2 and samtools
/// Splits each FASTQ into NUM_PARTS (default 4) parts, runs minimap2 concurrently, then merges BAMs
fn map_reads(config: &EvaluateConfig, outpath: &str, genome_size_bp: u64, log: &mut std::fs::File) -> Result<()> {
    use std::io::Write;

    // Pre-flight: verify tools are in PATH before touching any data files
    utils::require_tool("minimap2")
        .map_err(|e| anyhow::anyhow!("{}\n  Needed for: read-to-contig alignment", e))?;
    utils::require_tool("samtools")
        .map_err(|e| anyhow::anyhow!("{}\n  Needed for: BAM sorting and indexing", e))?;

    let preset = utils::get_minimap2_preset(&config.datatype);
    let valid_contig_fa = format!("{}valid_contig.fa", outpath);

    // Verify the contig FASTA written by Phase 2 is present
    utils::validate_input_file(&valid_contig_fa)
        .map_err(|e| anyhow::anyhow!("Phase 2 output missing — {}", e))?;

    const NUM_PARTS: usize = 4; // Split each FASTQ into 4 parts
    let mut all_bams = Vec::new();

    // Coverage estimation will be done after first split (no pre-scan)
    let mut total_input_reads = 0usize;
    let mut total_input_bases = 0u64;
    let mut keep_fraction: Option<f64> = None;
    let mut coverage_reported = false;

    for (read_idx, read_file) in config.read.iter().enumerate() {
        // Verify read file is accessible before starting alignment
        utils::validate_input_file(read_file)
            .map_err(|e| anyhow::anyhow!("Read file #{} — {}", read_idx + 1, e))?;

        info!("Splitting {} into {} parts", read_file, NUM_PARTS);

        // Create workspace for this read file's split parts
        let split_dir = format!("{}split_workspace_read_{}/", outpath, read_idx + 1);
        std::fs::create_dir_all(&split_dir)?;

        // Split the FASTQ into NUM_PARTS; per-part counts are returned during split (no re-read)
        let split_result = utils::split_fastq_gz(read_file, &split_dir, NUM_PARTS, keep_fraction)?;

        // On first read file, accumulate stats and estimate coverage
        if read_idx == 0 && !coverage_reported {
            total_input_reads += split_result.total_reads;
            total_input_bases += split_result.total_bases;

            // Estimate coverage from first file and decide keep fraction
            if genome_size_bp > 0 && total_input_bases > 0 {
                let estimated_coverage = total_input_bases as f64 / genome_size_bp as f64;
                if config.read_coverage > 0 && estimated_coverage > config.read_coverage as f64 {
                    keep_fraction = Some(config.read_coverage as f64 / estimated_coverage);
                    info!(
                        "Read file #{}: Estimated coverage = {:.2}x; subsampling to target {}x (keep fraction {:.4})",
                        read_idx + 1,
                        estimated_coverage,
                        config.read_coverage,
                        keep_fraction.unwrap(),
                    );
                    writeln!(
                        log,
                        "Read file #{}: Estimated coverage = {:.2}x; subsampling to target {}x (keep fraction {:.4})",
                        read_idx + 1,
                        estimated_coverage,
                        config.read_coverage,
                        keep_fraction.unwrap(),
                    )?;
                } else if total_input_bases > 0 {
                    info!(
                        "Read file #{}: Estimated coverage = {:.2}x; no subsampling needed (target {}x)",
                        read_idx + 1,
                        estimated_coverage,
                        config.read_coverage,
                    );
                    writeln!(
                        log,
                        "Read file #{}: Estimated coverage = {:.2}x; no subsampling needed (target {}x)",
                        read_idx + 1,
                        estimated_coverage,
                        config.read_coverage,
                    )?;
                }
            } else {
                warn!(
                    "Could not estimate read coverage (genome: {}, bases: {}); skipping subsampling",
                    genome_size_bp, total_input_bases,
                );
                writeln!(
                    log,
                    "Could not estimate read coverage (genome: {}, bases: {}); skipping subsampling",
                    genome_size_bp, total_input_bases,
                )?;
            }
            coverage_reported = true;

            // If keep_fraction is set, re-split with subsampling
            if keep_fraction.is_some() {
                std::fs::remove_dir_all(&split_dir)?;
                std::fs::create_dir_all(&split_dir)?;
                let split_result = utils::split_fastq_gz(read_file, &split_dir, NUM_PARTS, keep_fraction)?;
                info!(
                    "Subsampled {} from {} reads / {} bp to {} reads / {} bp",
                    read_file,
                    split_result.total_reads,
                    split_result.total_bases,
                    split_result.written_reads,
                    split_result.written_bases,
                );
                writeln!(
                    log,
                    "Subsampled {} from {} reads / {} bp to {} reads / {} bp",
                    read_file,
                    split_result.total_reads,
                    split_result.total_bases,
                    split_result.written_reads,
                    split_result.written_bases,
                )?;
                process_split_bams(&split_result, config, NUM_PARTS, &valid_contig_fa, &split_dir, read_idx, outpath, log)?;
                all_bams.push(format!("{}read_to_contig_{}.bam", outpath, read_idx + 1));
            } else {
                process_split_bams(&split_result, config, NUM_PARTS, &valid_contig_fa, &split_dir, read_idx, outpath, log)?;
                all_bams.push(format!("{}read_to_contig_{}.bam", outpath, read_idx + 1));
            }
        } else if read_idx == 0 {
            // Should not get here, but handle it
            process_split_bams(&split_result, config, NUM_PARTS, &valid_contig_fa, &split_dir, read_idx, outpath, log)?;
            all_bams.push(format!("{}read_to_contig_{}.bam", outpath, read_idx + 1));
        } else {
            // Subsequent read files: use same keep_fraction as first file
            info!(
                "Read file #{} with keep_fraction = {}",
                read_idx + 1,
                keep_fraction.map(|f| format!("{:.4}", f)).unwrap_or_else(|| "None".to_string()),
            );
            process_split_bams(&split_result, config, NUM_PARTS, &valid_contig_fa, &split_dir, read_idx, outpath, log)?;
            all_bams.push(format!("{}read_to_contig_{}.bam", outpath, read_idx + 1));
        }
    }

    // Helper to process split files and create merged BAM
    fn process_split_bams(
        split_result: &utils::SplitFastqResult,
        config: &EvaluateConfig,
        num_parts: usize,
        valid_contig_fa: &str,
        split_dir: &str,
        read_idx: usize,
        outpath: &str,
        log: &mut std::fs::File,
    ) -> Result<()> {
        use std::io::Write;

        info!("Reads per split part for read file {}:", read_idx + 1);
        writeln!(log, "Reads per split part for read file {}:", read_idx + 1)?;
        for (part_idx, &part_reads) in split_result.part_counts.iter().enumerate() {
            info!("  part_{:02}.fastq.gz: {} reads", part_idx, part_reads);
            writeln!(log, "  part_{:02}.fastq.gz: {} reads", part_idx, part_reads)?;
        }

        info!("Aligning {} parts (read file {}/{})", num_parts, read_idx + 1, config.read.len());

        // Process all parts in parallel using rayon
        let part_bams: Vec<String> = split_result.paths.par_iter()
            .enumerate()
            .map(|(part_idx, split_fastq)| {
                let bam_out = format!("{}part_{:02}.bam", split_dir, part_idx);

                // Run minimap2 + samtools for this part
                let minimap = std::process::Command::new("minimap2")
                    .arg("-a")
                    .arg("-Q")
                    .arg("-N")
                    .arg("1")
                    .arg("-I")
                    .arg("10G")
                    .arg("-t")
                    .arg((config.thread / num_parts).max(1).to_string())
                    .arg(valid_contig_fa)
                    .arg(split_fastq)
                    .stdout(std::process::Stdio::piped())
                    .spawn()
                    .map_err(|e| anyhow::anyhow!("Failed to start minimap2 on part {}: {}", part_idx, e))?;

                let samtools = std::process::Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg((config.thread / num_parts).max(1).to_string())
                    .arg("-o")
                    .arg(&bam_out)
                    .stdin(minimap.stdout.unwrap())
                    .output()
                    .map_err(|e| anyhow::anyhow!("Failed to start samtools sort on part {}: {}", part_idx, e))?;

                if !samtools.status.success() {
                    return Err(anyhow::anyhow!("samtools sort failed on part {}", part_idx));
                }

                Ok(bam_out)
            })
            .collect::<Result<Vec<_>>>()?;

        // Merge BAMs from all parts for this read file
        let read_merged_bam = format!("{}read_to_contig_{}.bam", outpath, read_idx + 1);

        if part_bams.len() > 1 {
            let merge = std::process::Command::new("samtools")
                .arg("merge")
                .arg("-o")
                .arg(&read_merged_bam)
                .args(&part_bams)
                .output()
                .map_err(|e| anyhow::anyhow!("Failed to merge BAM parts: {}", e))?;

            if !merge.status.success() {
                anyhow::bail!("samtools merge failed");
            }
            info!("Merged {} BAM parts to {}", part_bams.len(), read_merged_bam);
        } else if part_bams.len() == 1 {
            info!("Renaming single BAM from {} to {}", part_bams[0], read_merged_bam);
            std::fs::rename(&part_bams[0], &read_merged_bam)
                .map_err(|e| anyhow::anyhow!("Could not rename BAM from {} to {}: {}", part_bams[0], read_merged_bam, e))?;
        }

        // Clean up split workspace
        std::fs::remove_dir_all(split_dir)
            .map_err(|e| anyhow::anyhow!("Could not clean up split workspace: {}", e))?;

        Ok(())
    }

    // Merge BAMs from all read files (if multiple)
    if all_bams.len() > 1 {
        info!("Merging BAM files from {} read files", all_bams.len());
        let merge = std::process::Command::new("samtools")
            .arg("merge")
            .arg(format!("{}read_to_contig.bam", outpath))
            .args(&all_bams)
            .output()
            .map_err(|e| anyhow::anyhow!("Failed to merge read file BAMs: {}", e))?;

        if !merge.status.success() {
            let stderr = String::from_utf8_lossy(&merge.stderr);
            anyhow::bail!("samtools merge failed: {}", stderr);
        }

        // Cleanup individual read file BAMs
        for bam in all_bams {
            std::fs::remove_file(&bam)
                .map_err(|e| anyhow::anyhow!("Could not remove intermediate BAM '{}': {}", bam, e))?;
        }
    } else if all_bams.len() == 1 {
        // Only one read file, rename its BAM
        std::fs::rename(&all_bams[0], format!("{}read_to_contig.bam", outpath))
            .map_err(|e| anyhow::anyhow!("Could not rename final BAM: {}", e))?;
    }

    // Index the final BAM
    info!("Indexing BAM file");
    let index = std::process::Command::new("samtools")
        .arg("index")
        .arg(format!("{}read_to_contig.bam", outpath))
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to start samtools index: {}", e))?;

    if !index.status.success() {
        let stderr = String::from_utf8_lossy(&index.stderr);
        anyhow::bail!("samtools index failed: {}", stderr);
    }

    // Count mapped reads in the final BAM using samtools stats
    let stats_output = std::process::Command::new("samtools")
        .arg("stats")
        .arg(format!("{}read_to_contig.bam", outpath))
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to run samtools stats: {}", e))?;

    if stats_output.status.success() {
        let stats_str = String::from_utf8_lossy(&stats_output.stdout);
        // Parse SN (summary number) lines to get read counts
        let mut total_reads_in_bam = 0;
        let mut mapped_reads = 0;
        
        for line in stats_str.lines() {
            if line.starts_with("SN\tsequences:") {
                total_reads_in_bam = line.split('\t').last()
                    .and_then(|s| s.trim().parse::<u64>().ok())
                    .unwrap_or(0);
            } else if line.starts_with("SN\treads mapped:") {
                mapped_reads = line.split('\t').last()
                    .and_then(|s| s.trim().parse::<u64>().ok())
                    .unwrap_or(0);
            }
        }

        info!("Final BAM statistics: {} total reads, {} mapped", total_reads_in_bam, mapped_reads);
        writeln!(log, "Final BAM statistics: {} total reads, {} mapped", total_reads_in_bam, mapped_reads)?;
    }

    info!("Read mapping and BAM merge complete");
    Ok(())
}

/// Helper function: wraps public run() function for CLI compatibility
pub fn evaluate(config: EvaluateConfig) -> Result<()> {
    run(config)
}
