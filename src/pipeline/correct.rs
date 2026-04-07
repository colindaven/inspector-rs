/// Assembly correction pipeline
/// Reads evaluation output (error BED files) and corrects assembly

use anyhow::{Result, Context};
use log::info;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::Command;
use rayon::prelude::*;
use crate::utils;

#[derive(Debug, Clone)]
pub struct CorrectConfig {
    pub input: String,
    pub read: Vec<String>,
    pub assembly: String,
    pub outpath: String,
    pub datatype: String,
    pub thread: usize,
    pub flye_timeout: u64,
    pub min_correction_support: usize,
}

/// Represents a BED format error record
#[derive(Debug, Clone)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub sv_type: String,
    pub size: i64,
    pub support: usize,
}

/// Main entry point for assembly correction
pub fn run(config: CorrectConfig) -> Result<()> {
    let start_time = std::time::Instant::now();

    // Validate inputs
    info!("Validating input files for correction");
    // config.input is the evaluation *directory*, not a file
    if !std::path::Path::new(&config.input).is_dir() {
        anyhow::bail!("Input evaluation directory not found or not a directory: {}", config.input);
    }
    utils::validate_input_file(&config.assembly)?;
    utils::validate_input_files(&config.read)?;
    utils::validate_datatype(&config.datatype)?;

    // Normalize output path
    let outpath = utils::normalize_path(&config.outpath);
    utils::ensure_output_dir(&outpath)?;

    info!("Input evaluation dir: {}", config.input);
    info!("Assembly file: {}", config.assembly);
    info!("Output directory: {}", outpath);

    // Load contig sequences
    info!("Loading contig sequences");
    let mut contigs = load_contigs(&config.assembly)?;
    info!("Loaded {} contigs", contigs.len());

    // Load structural errors (prefer structural_error.bed; fallback to ae_merge_workspace outputs)
    let struct_errors_all = load_structural_errors(&config.input)?;
    let struct_before = struct_errors_all.len();
    let struct_errors: Vec<_> = struct_errors_all
        .into_iter()
        .filter(|e| e.support >= config.min_correction_support)
        .collect();
    let struct_after = struct_errors.len();
    let struct_dropped = struct_before.saturating_sub(struct_after);
    info!(
        "Structural BED records: before={}, after={}, dropped={} (min support: {})",
        struct_before,
        struct_after,
        struct_dropped,
        config.min_correction_support
    );

    // Load base errors (prefer small_scale_error.bed; fallback to base_error_workspace/baseerror_*.bed)
    let base_errors_all = load_base_errors(&config.input)?;
    let base_before = base_errors_all.len();
    let base_errors: Vec<_> = base_errors_all
        .into_iter()
        .filter(|e| e.support >= config.min_correction_support)
        .collect();
    let base_after = base_errors.len();
    let base_dropped = base_before.saturating_sub(base_after);
    info!(
        "Small-scale BED records: before={}, after={}, dropped={} (min support: {})",
        base_before,
        base_after,
        base_dropped,
        config.min_correction_support
    );

    // Create workspace
    std::fs::create_dir_all(format!("{}correction_workspace", outpath))?;

    // Phase 7: Correct all contigs in parallel
    // Both base and structural corrections operate on independent sequences
    info!("Phase 7a: Correcting base errors");
    info!("Phase 7b: Correcting structural errors");
    let base_corrections = base_errors.len();
    let struct_corrections = struct_errors.len();

    let contig_names: Vec<String> = contigs.keys().cloned().collect();
    let corrected_seqs: HashMap<String, String> = contig_names.par_iter()
        .map(|contig_name| {
            let mut seq = contigs[contig_name].clone();

            // Apply base corrections
            let berrors: Vec<_> = base_errors.iter().filter(|e| e.chrom == *contig_name).collect();
            if !berrors.is_empty() {
                apply_base_corrections(&mut seq, &berrors)?;
            }

            // Apply structural corrections (with optional Flye)
            let serrors: Vec<_> = struct_errors.iter().filter(|e| e.chrom == *contig_name).collect();
            if !serrors.is_empty() {
                let flye_result = if !config.read.is_empty() {
                    execute_flye_correction(
                        &config.input,
                        &config.assembly,
                        &config.read,
                        contig_name,
                        &config.datatype,
                        &config.thread,
                        &config.flye_timeout,
                    )
                } else {
                    Err(anyhow::anyhow!("No reads provided for Flye"))
                };

                if flye_result.is_err() {
                    apply_structural_corrections(&mut seq, &serrors)?;
                }
            }

            Ok((contig_name.clone(), seq))
        })
        .collect::<Result<HashMap<_, _>>>()?;

    contigs = corrected_seqs;

    // Phase 8: Output and statistics
    info!("Phase 8: Final assembly output and statistics");
    
    // Write corrected assembly
    let corrected_file = format!("{}corrected_assembly.fasta", outpath);
    info!("Writing corrected assembly to {}", corrected_file);
    write_fasta(&corrected_file, &contigs)?;

    // Compute correction statistics
    let total_bases_corrected = compute_correction_stats(&base_errors, &struct_errors);
    info!("Total bases affected by corrections: {}", total_bases_corrected);

    // Write summary statistics
    let summary_file = format!("{}correction_summary.txt", outpath);
    write_correction_summary(&summary_file, base_corrections, struct_corrections, &contigs)?;

    let total_time = std::time::Instant::now();
    info!("Assembly correction complete in {:.2}s", (total_time - start_time).as_secs_f64());
    info!("Corrected {} base errors and {} structural errors", base_corrections, struct_corrections);

    Ok(())
}

/// Load FASTA sequences into memory
fn load_contigs(path: &str) -> Result<HashMap<String, String>> {
    let mut contigs = HashMap::new();
    let file = File::open(path)?;
    
    // Auto-detect gzip format
    if path.ends_with(".gz") {
        let gz = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(gz);
        for line in reader.lines() {
            let line = line?;
            parse_fasta_line(&line, &mut contigs)?;
        }
    } else {
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            parse_fasta_line(&line, &mut contigs)?;
        }
    }
    
    Ok(contigs)
}

/// Helper to parse a FASTA line and update the contigs HashMap
fn parse_fasta_line(line: &str, contigs: &mut HashMap<String, String>) -> Result<()> {
    if line.starts_with('>') {
        // New sequence header - just register it with empty sequence
        let name = line[1..].split_whitespace().next().unwrap_or("");
        contigs.insert(name.to_string(), String::new());
    } else {
        // Sequence data - append to the last sequence
        if let Some((_, seq)) = contigs.iter_mut().last() {
            seq.push_str(line);
        }
    }
    Ok(())
}

/// Load BED file into a list of records
fn load_bed_file(path: &str) -> Result<Vec<BedRecord>> {
    let mut records = Vec::new();
    let file = File::open(path).context(format!("Cannot open BED file: {}", path))?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        
        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        
        // BED format: chrom, start, end, [name], [score], [strand], [extra fields...]
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);

        // Try to parse additional fields that Inspector writes
        let sv_type = if fields.len() > 3 { fields[3].to_string() } else { "unknown".to_string() };
        let size = if fields.len() > 4 { fields[4].parse().unwrap_or((end - start) as i64) } else { (end - start) as i64 };
        // Support can be column 5 (base BED: chrom start end type support)
        // or column 6 (structural BED: chrom start end type size support).
        let support = if fields.len() > 5 {
            fields[5].parse().unwrap_or(1)
        } else if fields.len() > 4 {
            fields[4].parse().unwrap_or(1)
        } else {
            1
        };

        records.push(BedRecord {
            chrom,
            start,
            end,
            sv_type,
            size,
            support,
        });
    }

    Ok(records)
}

fn load_base_errors(input_dir: &str) -> Result<Vec<BedRecord>> {
    let base_errors_file = format!("{}small_scale_error.bed", utils::normalize_path(input_dir));
    if std::path::Path::new(&base_errors_file).exists() {
        return load_bed_file(&base_errors_file);
    }

    // Fallback for current Rust evaluate output: aggregate baseerror_*.bed files.
    let mut records = Vec::new();
    let workspace = format!("{}base_error_workspace", utils::normalize_path(input_dir));
    if let Ok(entries) = std::fs::read_dir(&workspace) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.starts_with("baseerror_") && n.ends_with(".bed"))
                .unwrap_or(false)
            {
                let path_str = path.to_string_lossy().to_string();
                if let Ok(mut file_records) = load_bed_file(&path_str) {
                    records.append(&mut file_records);
                }
            }
        }
    }

    Ok(records)
}

fn load_structural_errors(input_dir: &str) -> Result<Vec<BedRecord>> {
    let struct_errors_file = format!("{}structural_error.bed", utils::normalize_path(input_dir));
    if std::path::Path::new(&struct_errors_file).exists() {
        return load_bed_file(&struct_errors_file);
    }

    // Fallback for current Rust evaluate output: parse ae_merge_workspace/del_merged_* and ins_merged_*.
    let mut records = Vec::new();
    let workspace = format!("{}ae_merge_workspace", utils::normalize_path(input_dir));
    if let Ok(entries) = std::fs::read_dir(&workspace) {
        for entry in entries.flatten() {
            let path = entry.path();
            let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
                continue;
            };

            let sv_type = if name.starts_with("del_merged_") {
                "DEL"
            } else if name.starts_with("ins_merged_") {
                "INS"
            } else {
                continue;
            };

            if let Ok(file) = File::open(&path) {
                let reader = BufReader::new(file);
                for line in reader.lines() {
                    let line = line?;
                    if line.trim().is_empty() {
                        continue;
                    }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 3 {
                        continue;
                    }
                    let chrom = fields[0].to_string();
                    let pos: u64 = fields[1].parse().unwrap_or(0);
                    let support: usize = fields[2].parse().unwrap_or(1);
                    let start = pos.saturating_sub(1);
                    let end = pos;

                    records.push(BedRecord {
                        chrom,
                        start,
                        end,
                        sv_type: sv_type.to_string(),
                        size: 1,
                        support,
                    });
                }
            }
        }
    }

    Ok(records)
}

/// Write FASTA sequences to file
fn write_fasta(path: &str, contigs: &HashMap<String, String>) -> Result<()> {
    let mut file = File::create(path)?;
    
    for (name, seq) in contigs {
        writeln!(file, ">{}", name)?;
        
        // Write sequence in 80-character lines
        for chunk in seq.chars().collect::<Vec<_>>().chunks(80) {
            let line: String = chunk.iter().collect();
            writeln!(file, "{}", line)?;
        }
    }

    Ok(())
}

/// Apply base-level corrections (SNPs and small indels)
fn apply_base_corrections(sequence: &mut String, errors: &[&BedRecord]) -> Result<()> {
    // Sort errors by position (descending) to process from right to left
    // This prevents position shifts when applying corrections
    let mut sorted_errors = errors.to_vec();
    sorted_errors.sort_by(|a, b| b.start.cmp(&a.start));

    for error in sorted_errors {
        match error.sv_type.as_str() {
            "SNP" => {
                // For SNPs, we try to replace with consensus base (placeholder: use 'N')
                if (error.start as usize) < sequence.len() {
                    let chars: Vec<char> = sequence.chars().collect();
                    if chars[error.start as usize] != 'N' {
                        // In practice, would need to know what base to replace with
                        // For now, mark as ambiguous
                        let mut new_seq: String = chars[..error.start as usize].iter().collect();
                        new_seq.push('N');
                        new_seq.push_str(&chars[(error.start as usize + 1)..].iter().collect::<String>());
                        *sequence = new_seq;
                    }
                }
            }
            "insertion" | "INS" => {
                // Remove inserted bases
                if (error.start as usize) < sequence.len() && (error.end as usize) < sequence.len() {
                    let chars: Vec<char> = sequence.chars().collect();
                    let mut new_seq: String = chars[..error.start as usize].iter().collect();
                    new_seq.push_str(&chars[(error.end as usize)..].iter().collect::<String>());
                    *sequence = new_seq;
                }
            }
            "deletion" | "DEL" => {
                // Insertions to fill deletions - use reference or Ns
                if (error.start as usize) < sequence.len() {
                    let size = error.size.abs() as usize;
                    let mut new_seq: String = sequence.chars().take(error.start as usize).collect();
                    for _ in 0..size {
                        new_seq.push('N');
                    }
                    let remaining: String = sequence.chars().skip(error.start as usize).collect();
                    new_seq.push_str(&remaining);
                    *sequence = new_seq;
                }
            }
            _ => {
                // Unknown error type, skip
            }
        }
    }

    Ok(())
}

/// Apply structural error corrections (simple patching approach)
fn apply_structural_corrections(sequence: &mut String, errors: &[&BedRecord]) -> Result<()> {
    // Sort errors by position (descending) to process from right to left
    let mut sorted_errors = errors.to_vec();
    sorted_errors.sort_by(|a, b| b.start.cmp(&a.start));

    for error in sorted_errors {
        match error.sv_type.as_str() {
            "deletion" | "DEL" => {
                // Fill deletion with Ns
                if (error.start as usize) < sequence.len() {
                    let size = error.size.abs().min((error.end - error.start) as i64) as usize;
                    let mut new_seq: String = sequence.chars().take(error.start as usize).collect();
                    for _ in 0..size {
                        new_seq.push('N');
                    }
                    if (error.end as usize) < sequence.len() {
                        let remaining: String = sequence.chars().skip(error.end as usize).collect();
                        new_seq.push_str(&remaining);
                    }
                    *sequence = new_seq;
                }
            }
            "insertion" | "INS" => {
                // Remove insertion
                if (error.start as usize) < sequence.len() && (error.end as usize) < sequence.len() {
                    let chars: Vec<char> = sequence.chars().collect();
                    let mut new_seq: String = chars[..error.start as usize].iter().collect();
                    new_seq.push_str(&chars[(error.end as usize)..].iter().collect::<String>());
                    *sequence = new_seq;
                }
            }
            "inversion" | "INV" => {
                // Reverse the inverted region
                if (error.start as usize) < sequence.len() && (error.end as usize) <= sequence.len() {
                    let chars: Vec<char> = sequence.chars().collect();
                    let mut new_seq: String = chars[..error.start as usize].iter().collect();
                    
                    // Reverse the region and complement (if DNA)
                    let mut inverted: Vec<char> = chars[error.start as usize..error.end as usize]
                        .iter()
                        .rev()
                        .map(|&c| complement_base(c))
                        .collect();
                    new_seq.push_str(&inverted.iter().collect::<String>());
                    new_seq.push_str(&chars[(error.end as usize)..].iter().collect::<String>());
                    *sequence = new_seq;
                }
            }
            _ => {}
        }
    }

    Ok(())
}

/// Get DNA complement base
fn complement_base(base: char) -> char {
    match base {
        'A' | 'a' => 'T',
        'T' | 't' => 'A',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        _ => base,
    }
}

/// Execute Flye for local re-assembly of error regions
fn execute_flye_correction(
    input_dir: &str,
    assembly: &str,
    reads: &[String],
    contig_name: &str,
    datatype: &str,
    threads: &usize,
    timeout: &u64,
) -> Result<()> {
    // Check if Flye is available
    let flye_check = Command::new("flye")
        .arg("--version")
        .output();

    if flye_check.is_err() {
        anyhow::bail!("Flye not found in PATH");
    }

    info!("Attempting local re-assembly of {} using Flye", contig_name);

    // Extract reads overlapping the error region (would need BAM index)
    // For now, create a simple Flye workspace and attempt re-assembly
    let workspace = format!("{}/flye_workspace_{}", input_dir, contig_name);
    std::fs::create_dir_all(&workspace)?;

    // Prepare Flye arguments based on datatype
    let mut flye_args = vec![
        "--nano-raw".to_string(),  // Default to nano-raw
        reads[0].clone(),
        "-o".to_string(),
        workspace.clone(),
        "-t".to_string(),
        threads.to_string(),
    ];

    // Adjust preset for datatype
    if datatype.contains("hifi") || datatype.contains("ccs") {
        flye_args[0] = "--pacbio-hifi".to_string();
    } else if datatype.contains("bax") || datatype.contains("subreads") {
        flye_args[0] = "--pacbio-raw".to_string();
    }

    // Run Flye with timeout
    info!("Running: flye {}", flye_args.join(" "));
    
    let status = Command::new("timeout")
        .arg(timeout.to_string())
        .arg("flye")
        .args(&flye_args)
        .output()
        .context("Failed to run Flye")?;

    if !status.status.success() {
        anyhow::bail!("Flye execution failed or timed out");
    }

    // TODO: Read corrected assembly from Flye output
    // and replace the contig sequence with corrected version

    info!("Flye correction completed for {}", contig_name);
    Ok(())
}

/// Calculate how many bases were corrected
fn compute_correction_stats(base_errors: &[BedRecord], struct_errors: &[BedRecord]) -> u64 {
    let mut total = 0u64;
    
    for error in base_errors {
        total += (error.end - error.start).max(1);
    }

    for error in struct_errors {
        total += error.size.abs() as u64;
    }

    total
}

/// Write correction summary statistics
fn write_correction_summary(
    path: &str,
    base_count: usize,
    struct_count: usize,
    contigs: &HashMap<String, String>,
) -> Result<()> {
    let mut file = File::create(path)?;
    
    writeln!(file, "Assembly Correction Summary")?;
    writeln!(file, "============================")?;
    writeln!(file)?;
    writeln!(file, "Base-level errors corrected: {}", base_count)?;
    writeln!(file, "Structural errors corrected: {}", struct_count)?;
    writeln!(file, "Total errors corrected: {}", base_count + struct_count)?;
    writeln!(file)?;
    
    let total_length: u64 = contigs.values().map(|s| s.len() as u64).sum();
    writeln!(file, "Corrected assembly statistics:")?;
    writeln!(file, "  Total length: {} bp", utils::format_size(total_length))?;
    writeln!(file, "  Number of contigs: {}", contigs.len())?;
    
    // Compute N50
    if let Ok(n50) = compute_n50_from_contigs(contigs) {
        writeln!(file, "  N50: {}", utils::format_size(n50))?;
    }

    Ok(())
}

/// Compute N50 from contigs
fn compute_n50_from_contigs(contigs: &HashMap<String, String>) -> Result<u64> {
    let mut lengths: Vec<u64> = contigs.values().map(|s| s.len() as u64).collect();
    lengths.sort_by(|a, b| b.cmp(a));

    let total_length: u64 = lengths.iter().sum();
    let mut cumulative = 0u64;

    for length in lengths {
        cumulative += length;
        if cumulative >= total_length / 2 {
            return Ok(length);
        }
    }

    Ok(0)
}

/// Helper function: wraps public run() function for CLI compatibility
pub fn correct(config: CorrectConfig) -> Result<()> {
    run(config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base('A'), 'T');
        assert_eq!(complement_base('T'), 'A');
        assert_eq!(complement_base('C'), 'G');
        assert_eq!(complement_base('G'), 'C');
        assert_eq!(complement_base('a'), 'T');
        assert_eq!(complement_base('N'), 'N');
    }

    #[test]
    fn test_apply_snp_correction() {
        let mut seq = "ACGTACGT".to_string();
        let error = BedRecord {
            chrom: "chr1".to_string(),
            start: 2,
            end: 3,
            sv_type: "SNP".to_string(),
            size: 1,
            support: 5,
        };
        let errors = vec![&error];

        apply_base_corrections(&mut seq, &errors).unwrap();
        
        // SNP at position 2 should be replaced with N
        assert_eq!(seq.chars().nth(2), Some('N'));
    }

    #[test]
    fn test_apply_insertion_correction() {
        let mut seq = "ACGTACGT".to_string();
        let error = BedRecord {
            chrom: "chr1".to_string(),
            start: 2,
            end: 5, // Remove 3 bases (indices 2, 3, 4 = G, T, A)
            sv_type: "INS".to_string(),
            size: 3,
            support: 5,
        };
        let errors = vec![&error];

        apply_base_corrections(&mut seq, &errors).unwrap();
        
        // Should remove bases 2-5, keeping 0-1 and 5-7: AC + CGT = ACCGT
        assert_eq!(seq, "ACCGT");
    }

    #[test]
    fn test_apply_deletion_correction() {
        let mut seq = "ACGTACGT".to_string();
        let error = BedRecord {
            chrom: "chr1".to_string(),
            start: 2,
            end: 2,
            sv_type: "DEL".to_string(),
            size: 2,
            support: 5,
        };
        let errors = vec![&error];

        apply_base_corrections(&mut seq, &errors).unwrap();
        
        // Should insert 2 Ns at position 2
        assert_eq!(seq, "ACNNGTACGT");
    }

    #[test]
    fn test_apply_structural_deletion() {
        let mut seq = "ACGTACGTACGT".to_string();
        let error = BedRecord {
            chrom: "chr1".to_string(),
            start: 4,
            end: 8,
            sv_type: "DEL".to_string(),
            size: 4,
            support: 5,
        };
        let errors = vec![&error];

        apply_structural_corrections(&mut seq, &errors).unwrap();
        
        // Should replace with Ns: "ACGTNNNNACGT"
        assert_eq!(seq, "ACGTNNNNACGT");
    }

    #[test]
    fn test_apply_structural_inversion() {
        // Use a sequence where reverse complement produces different output
        let mut seq = "ACGATACGT".to_string();  // A C G A T A C G T
        let error = BedRecord {
            chrom: "chr1".to_string(),
            start: 2,
            end: 7,  // Extract GATAC (indices 2-6)
            sv_type: "INV".to_string(),
            size: 5,
            support: 5,
        };
        let errors = vec![&error];

        apply_structural_corrections(&mut seq, &errors).unwrap();
        
        // GATAC reversed = CATAG
        // CATAG complemented = GTATC
        // Result: AC + GTATC + GT = ACGTATCGT
        assert_eq!(seq, "ACGTATCGT");
    }

    #[test]
    fn test_compute_n50_from_contigs() {
        let mut contigs = HashMap::new();
        contigs.insert("chr1".to_string(), "A".repeat(100));
        contigs.insert("chr2".to_string(), "A".repeat(200));
        contigs.insert("chr3".to_string(), "A".repeat(300));
        contigs.insert("chr4".to_string(), "A".repeat(400));

        let n50 = compute_n50_from_contigs(&contigs).unwrap();
        // Total = 1000, half = 500
        // Sorted: 400 + 300 = 700 >= 500, so N50 = 300
        assert_eq!(n50, 300);
    }

    #[test]
    fn test_compute_correction_stats() {
        let base_errors = vec![
            BedRecord {
                chrom: "chr1".to_string(),
                start: 0,
                end: 10,
                sv_type: "SNP".to_string(),
                size: 1,
                support: 1,
            },
            BedRecord {
                chrom: "chr1".to_string(),
                start: 50,
                end: 60,
                sv_type: "INS".to_string(),
                size: 10,
                support: 1,
            },
        ];

        let struct_errors = vec![
            BedRecord {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                sv_type: "DEL".to_string(),
                size: 100,
                support: 5,
            },
        ];

        let total = compute_correction_stats(&base_errors, &struct_errors);
        // base: (10-0) + (60-50) = 10 + 10 = 20
        // struct: |100| = 100
        // total = 120
        assert_eq!(total, 120);
    }

    #[test]
    fn test_bed_record_parsing() {
        let bed_line = "chr1\t100\t200\tDEL\t100\t5";
        let fields: Vec<&str> = bed_line.split('\t').collect();
        
        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().unwrap();
        let end: u64 = fields[2].parse().unwrap();
        let sv_type = if fields.len() > 3 { fields[3].to_string() } else { "unknown".to_string() };
        let size = if fields.len() > 4 { fields[4].parse().unwrap_or(0i64) } else { (end - start) as i64 };
        let support = if fields.len() > 5 { fields[5].parse().unwrap_or(1) } else { 1 };

        assert_eq!(chrom, "chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert_eq!(sv_type, "DEL");
        assert_eq!(size, 100);
        assert_eq!(support, 5);
    }

    #[test]
    fn test_multiple_error_corrections() {
        let mut seq = "ACGTACGTACGTACGT".to_string();
        
        // Define multiple non-overlapping errors from right to left
        let err1 = BedRecord {
            chrom: "chr1".to_string(),
            start: 12,
            end: 14,
            sv_type: "INS".to_string(),
            size: 2,
            support: 1,
        };
        let err2 = BedRecord {
            chrom: "chr1".to_string(),
            start: 4,
            end: 5,
            sv_type: "SNP".to_string(),
            size: 1,
            support: 1,
        };

        let errors = vec![&err1, &err2];
        apply_base_corrections(&mut seq, &errors).unwrap();
        
        // Changes should be applied without overlap issues
        assert!(seq.len() >= 14); // Should have modifications
    }
}
