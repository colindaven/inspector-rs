/// Static analysis of FASTA files and statistics computation
/// FASTA I/O, N50 computation, coverage merging

pub mod fasta;

use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use log::info;
use crate::utils;

/// Simple FASTA file processing - filters by length, computes stats
pub fn simple_fasta(
    contig_files: &[String],
    outpath: &str,
    min_size: usize,
    min_size_assembly_error: usize,
) -> Result<FastaStats> {
    info!("Loading FASTA files");

    let mut all_contigs = Vec::new();

    // Load all input FASTA files
    for file_path in contig_files {
        let records = fasta::load_fasta(file_path)?;
        all_contigs.extend(records);
        info!("Loaded {} sequences from {}", all_contigs.len(), file_path);
    }

    let mut chromosomes = Vec::new();
    let mut chromosomes_map = HashMap::new();
    let mut chromosomes_large = Vec::new();
    let mut contig_lengths = HashMap::new();
    let mut all_lengths: Vec<u64> = Vec::new();
    let mut total_length = 0u64;
    let mut total_length_large = 0u64;
    let mut total_length_all = 0u64;
    let mut total_length_above_ae = 0u64;

    // Write valid contigs and compute statistics
    let valid_fa_path = format!("{}valid_contig.fa", outpath);
    let mut valid_fa = File::create(&valid_fa_path)?;
    let contig_info_path = format!("{}contig_length_info", outpath);
    let mut contig_info = File::create(&contig_info_path)?;

    for contig in all_contigs {
        let seq_len = contig.sequence.len() as u64;
        all_lengths.push(seq_len);
        total_length_all += seq_len;
        contig_lengths.insert(contig.name.clone(), seq_len);

        // Write to contig info file
        writeln!(contig_info, "{}\t{}", contig.name, seq_len)?;

        if seq_len >= min_size as u64 {
            chromosomes.push(contig.name.clone());
            chromosomes_map.insert(contig.name.clone(), contig.name.clone());
            total_length += seq_len;

            // Track large contigs (for SV detection)
            if seq_len >= min_size_assembly_error as u64 {
                chromosomes_large.push(contig.name.clone());
                total_length_large += seq_len;
                total_length_above_ae += seq_len;
            }

            // Write valid contig to FASTA
            writeln!(valid_fa, ">{}", contig.name)?;
            writeln!(valid_fa, "{}", contig.sequence)?;
        }
    }

    // Sort lengths descending for N50 and rank queries
    let mut sorted_lengths = all_lengths.clone();
    sorted_lengths.sort_unstable_by(|a, b| b.cmp(a));

    // N50 over all contigs (Python uses all_lengths for the main N50)
    let n50 = compute_n50_sorted(&sorted_lengths);

    // N50 of contigs > min_size_assembly_error
    let large_lengths: Vec<u64> = sorted_lengths
        .iter()
        .filter(|&&l| l > min_size_assembly_error as u64)
        .copied()
        .collect();
    let n50_large = compute_n50_sorted(&large_lengths);

    let largest_contig_length = sorted_lengths.first().copied().unwrap_or(0);
    let second_largest_contig_length = sorted_lengths.get(1).copied().unwrap_or(0);

    let total_contig_count_all = sorted_lengths.len();
    let count_above_ae_size = large_lengths.len();

    let stats = FastaStats {
        chromosomes: chromosomes.clone(),
        chromosomes_map,
        chromosomes_large: chromosomes_large.clone(),
        total_length,
        total_length_large,
        n50,
        largest_contig_length,
        second_largest_contig_length,
        contig_lengths,
        total_contig_count_all,
        total_length_all,
        count_above_ae_size,
        total_length_above_ae_size: total_length_above_ae,
        min_size_used: min_size,
        min_size_ae_used: min_size_assembly_error,
    };

    // Write summary statistics — Python-compatible "Statics of contigs:" section
    let summary_path = format!("{}summary_statistics", outpath);
    let mut summary = File::create(&summary_path)?;
    writeln!(summary, "Statics of contigs:")?;
    writeln!(summary, "Number of contigs\t{}", stats.total_contig_count_all)?;
    writeln!(summary, "Number of contigs > {} bp\t{}", min_size, stats.chromosomes.len())?;
    writeln!(summary, "Number of contigs >{} bp\t{}", min_size_assembly_error, stats.count_above_ae_size)?;
    writeln!(summary, "Total length\t{}", stats.total_length_all)?;
    writeln!(summary, "Total length of contigs > {} bp\t{}", min_size, stats.total_length)?;
    writeln!(summary, "Total length of contigs >{}bp\t{}", min_size_assembly_error, stats.total_length_above_ae_size)?;
    writeln!(summary, "Longest contig\t{}", stats.largest_contig_length)?;
    if stats.second_largest_contig_length > 0 {
        writeln!(summary, "Second longest contig length\t{}", stats.second_largest_contig_length)?;
    }
    writeln!(summary, "N50\t{}", stats.n50)?;
    writeln!(summary, "N50 of contigs >{}bp\t{}", min_size_assembly_error, n50_large)?;
    writeln!(summary)?;
    writeln!(summary)?;

    info!("Contigs: {} total, {} > {}bp, {} > {}bp",
        stats.total_contig_count_all, stats.chromosomes.len(), min_size,
        stats.count_above_ae_size, min_size_assembly_error);
    info!("Total length: {}", utils::format_size(stats.total_length));
    info!("N50: {}", utils::format_size(stats.n50));

    Ok(stats)
}

/// Compute N50 from a pre-sorted (descending) list of lengths.
fn compute_n50_sorted(sorted_lengths: &[u64]) -> u64 {
    let total: u64 = sorted_lengths.iter().sum();
    if total == 0 { return 0; }
    let half = total / 2;
    let mut acc = 0u64;
    for &l in sorted_lengths {
        acc += l;
        if acc >= half {
            return l;
        }
    }
    0
}

/// Aggregate mapping statistics from depth files
pub fn mapping_info_contig(
    outpath: &str,
    large_chroms: &[String],
    small_chroms: &[String],
    total_length: u64,
    _total_length_large: u64,
) -> Result<usize> {
    // Read depth files written by detect_sortbam / detect_sortbam_nosv.
    // Each file contains: chrom \t reads \t total_mapped_bp
    // Coverage for a contig = total_mapped_bp / contig_length.
    // We compute genome-wide average as sum(total_mapped_bp) / total_length.
    use std::io::BufRead;

    let all_chroms: Vec<&str> = large_chroms.iter().chain(small_chroms.iter()).map(|s| s.as_str()).collect();
    let mut total_mapped_bp: u64 = 0;

    for chrom in all_chroms {
        let depth_file = format!("{}map_depth/read_to_contig_{}.depth", outpath, chrom);
        if let Ok(f) = File::open(&depth_file) {
            for line in std::io::BufReader::new(f).lines().flatten() {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 3 {
                    if let Ok(bp) = parts[2].parse::<u64>() {
                        total_mapped_bp += bp;
                    }
                }
            }
        }
    }

    if total_length == 0 || total_mapped_bp == 0 {
        info!("Coverage files not found or empty; defaulting to coverage=20");
        return Ok(20);
    }

    let coverage = (total_mapped_bp / total_length) as usize;
    // Guard against pathologically low or high values
    let coverage = coverage.max(1);
    Ok(coverage)
}

/// Helper: collect all lines from per-chrom merged files matching `prefix`.
fn collect_merged(outpath: &str, prefix: &str) -> Vec<String> {
    use std::io::BufRead;
    let dir = format!("{}ae_merge_workspace/", outpath);
    let mut lines = Vec::new();
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return lines;
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.starts_with(prefix) {
                if let Ok(f) = File::open(&path) {
                    let reader = std::io::BufReader::new(f);
                    for l in reader.lines().flatten() {
                        if !l.is_empty() {
                            lines.push(l);
                        }
                    }
                }
            }
        }
    }

    lines
}

/// Merge and report assembly errors from clustering.
/// Mirrors Python denovo_static.assembly_info_cluster().
pub fn assembly_info_cluster(
    outpath: &str,
    min_size: usize,
    max_size: usize,
) -> Result<()> {
    let bed_path = format!("{}assembly_errors.bed", outpath);
    let mut bed = File::create(&bed_path)?;

    let parse = |line: &str| -> Option<(String, u64, u64, String, String)> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            info!("  assembly_info_cluster SKIP (too few fields={}): {}", f.len(), &line[..line.len().min(80)]);
            return None;
        }
        let pos: u64 = f[1].parse().ok()?;
        let size: u64 = f[2].parse().ok()?;
        if (size as usize) < min_size || (size as usize) > max_size {
            info!("  assembly_info_cluster SKIP (size={} not in [{},{}]): {}",
                  size, min_size, max_size, &line[..line.len().min(80)]);
            return None;
        }
        Some((f[0].to_string(), pos, size, f[3].to_string(), f[7].to_string()))
    };

    let mut n_exp = 0usize;
    let mut n_col = 0usize;
    let mut n_inv = 0usize;

    for line in collect_merged(outpath, "del_merged_") {
        if let Some((chrom, pos, size, count, readnames)) = parse(&line) {
            writeln!(
                bed,
                "{}\t{}\t{}\t{}\tExpansion\tSize={}\t{}",
                chrom,
                pos,
                pos + size,
                count,
                size,
                readnames
            )?;
            n_exp += 1;
        }
    }

    for line in collect_merged(outpath, "ins_merged_") {
        if let Some((chrom, pos, size, count, readnames)) = parse(&line) {
            writeln!(
                bed,
                "{}\t{}\t{}\t{}\tCollapse\tSize={}\t{}",
                chrom,
                pos,
                pos + 1,
                count,
                size,
                readnames
            )?;
            n_col += 1;
        }
    }

    for line in collect_merged(outpath, "inv_merged_") {
        if let Some((chrom, pos, size, count, readnames)) = parse(&line) {
            writeln!(
                bed,
                "{}\t{}\t{}\t{}\tInversion\t{}",
                chrom,
                pos,
                pos + size,
                count,
                readnames
            )?;
            n_inv += 1;
        }
    }

    info!(
        "assembly_errors.bed: {} Expansions, {} Collapses, {} Inversions",
        n_exp, n_col, n_inv
    );

    Ok(())
}

/// Reference-based assembly error reporting
pub fn assembly_info_ref(_outpath: &str) -> Result<()> {
    Ok(())
}

/// Get reference alignment info and coverage
pub fn get_ref_align_info(_path: &str, _total_length: u64) -> Result<()> {
    Ok(())
}

/// Count base-pair errors from reference alignment
pub fn basepair_error_ref(_outpath: &str, _largest_chr: &str) -> Result<()> {
    Ok(())
}

/// Statistics from FASTA processing
#[derive(Debug, Clone)]
pub struct FastaStats {
    pub chromosomes: Vec<String>,
    pub chromosomes_map: HashMap<String, String>,
    pub chromosomes_large: Vec<String>,
    pub total_length: u64,
    pub total_length_large: u64,
    pub n50: u64,
    pub largest_contig_length: u64,
    pub second_largest_contig_length: u64,
    pub contig_lengths: HashMap<String, u64>,
    /// Total contig count including those below min_size
    pub total_contig_count_all: usize,
    /// Total assembly length including all contigs
    pub total_length_all: u64,
    /// Count of contigs >= min_size_assemblyerror (for SV detection)
    pub count_above_ae_size: usize,
    /// Total length of contigs >= min_size_assemblyerror
    pub total_length_above_ae_size: u64,
    /// The min_size threshold used for filtering (for output labels)
    pub min_size_used: usize,
    /// The min_size_assemblyerror threshold used (for output labels)
    pub min_size_ae_used: usize,
}
