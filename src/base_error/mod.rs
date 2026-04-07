/// Base-level error detection (SNPs and small indels)
/// Parses pileup and detects errors using statistical tests

use std::process::{Command, Stdio};
use std::io::{BufRead, BufReader};
use std::fs::File;
use anyhow::{Result, Context};
use log::{debug, info};
use regex::Regex;

/// Detect SNVs and small indels from pileup
pub fn get_snv(
    path: &str,
    chrom: &str,
    mincount: usize,
    maxcov: usize,
    mindepth: Option<usize>,
) -> Result<()> {
    info!("Detecting SNVs/indels for {}", chrom);

    let bam_file = format!("{}read_to_contig.bam", path);
    let fa_file = format!("{}valid_contig.fa", path);
    let workspace = format!("{}base_error_workspace/", path);
    let output_file = format!("{}baseerror_{}.bed", workspace, chrom);

    // Run samtools mpileup
    let mpileup = Command::new("samtools")
        .arg("mpileup")
        .arg("-Q")
        .arg("0")
        .arg(&bam_file)
        .arg("-r")
        .arg(chrom)
        .arg("-o")
        .arg(format!("{}base_{}.pileup", workspace, chrom))
        .arg("-f")
        .arg(&fa_file)
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to run samtools mpileup: {}\n  BAM: {}\n  Reference: {}\n  Contig: {}\n  Hint: ensure samtools is in PATH and both files exist", e, bam_file, fa_file, chrom))?;

    if !mpileup.status.success() {
        anyhow::bail!("samtools mpileup failed for contig '{}':\n  BAM: {}\n  Reference: {}\n  stderr: {}",
            chrom, bam_file, fa_file, String::from_utf8_lossy(&mpileup.stderr));
    }

    // Parse pileup file
    let pileup_file = format!("{}base_{}.pileup", workspace, chrom);
    let file = File::open(&pileup_file)?;
    let reader = BufReader::new(file);

    let mut output = File::create(&output_file)?;
    use std::io::Write;

    let mindepth_val = mindepth.unwrap_or(maxcov / 10);
    let mut valid_bases = 0u64;
    let mut error_count = 0u64;

    // Regex for parsing pileup (simplified)
    let indent_re = Regex::new(r"\+\d+").ok();
    let del_re = Regex::new(r"-\d+").ok();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() < 5 {
            continue;
        }

        let ref_base = fields[2];
        let depth_str = fields[3];
        let bases_str = fields[4];

        let Ok(depth) = depth_str.parse::<usize>() else {
            continue;
        };

        if depth < mincount || depth > maxcov {
            continue;
        }

        // Count valid bases
        if ref_base != "N" && depth >= mindepth_val {
            valid_bases += 1;
        }

        // Parse base string: . or , = match, ACGT = mismatch, * = deletion
        let mut bases = bases_str.to_string();
        bases = bases.replace(",", ".");
        bases = bases.replace("a", "A");
        bases = bases.replace("t", "T");
        bases = bases.replace("c", "C");
        bases = bases.replace("g", "G");

        // Remove read start/end markers (^. pattern)
        if let Ok(re) = Regex::new(r"\^.") {
            bases = re.replace_all(&bases, "").to_string();
        }

        let deletions = bases.matches('*').count();
        let effective_depth = depth - deletions;

        if effective_depth == 0 {
            continue;
        }

        let min_support = ((depth as f64 * 0.2).max(mincount as f64)) as usize;

        // Detect insertions
        let insertions = if let Some(ref re) = indent_re {
            re.find_iter(&bases).count()
        } else {
            0
        };

        let pos: usize = fields[1].parse().unwrap_or(0);
        let bed_start = pos.saturating_sub(1); // pileup is 1-based; BED is 0-based

        if insertions >= min_support {
            debug!("Found insertion at {}", fields[1]);
            error_count += 1;
            writeln!(output, "{}\t{}\t{}\tinsertion\t{}", fields[0], bed_start, pos, insertions)?;
        }

        // Detect deletions (via * characters)
        if deletions >= min_support {
            debug!("Found deletion at {}", fields[1]);
            error_count += 1;
            writeln!(output, "{}\t{}\t{}\tdeletion\t{}", fields[0], bed_start, pos, deletions)?;
        }

        // Detect SNPs (mismatched bases)
        let mut base_counts = std::collections::HashMap::new();
        for c in bases.chars() {
            if c != '.' && c != '-' && c != '+' && c != '*' && c != '^' {
                *base_counts.entry(c).or_insert(0) += 1;
            }
        }

        for (_alt_base, count) in base_counts {
            if count >= min_support && count < effective_depth {
                // This is a SNP
                let freq = count as f64 / effective_depth as f64;
                if freq > 0.1 && freq < 0.9 { // Between 10% and 90% frequency
                    debug!("Found SNP at {}", fields[1]);
                    error_count += 1;
                    writeln!(output, "{}\t{}\t{}\tSNP\t{}", fields[0], bed_start, pos, count)?;
                }
            }
        }
    }

    // Write valid base count
    let validbase_file = format!("{}validbase", workspace);
    let mut validbase = File::create(&validbase_file)?;
    writeln!(validbase, "{}", valid_bases)?;

    info!("Found {} errors in {}", error_count, chrom);
    Ok(())
}

/// Count and filter base errors
pub fn count_base_errors(
    path: &str,
    ctg_total_length: u64,
    _datatype: &str,
    _ave_depth: usize,
) -> Result<usize> {
    debug!("Counting base errors");

    // Aggregate errors from all chromsomes
    let mut total_errors = 0usize;

    // List all baseerror_*.bed files in workspace
    let workspace = format!("{}base_error_workspace/", path);
    if let Ok(entries) = std::fs::read_dir(&workspace) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.starts_with("baseerror_") && n.ends_with(".bed"))
                .unwrap_or(false)
            {
                if let Ok(file) = File::open(&path) {
                    let reader = BufReader::new(file);
                    for line in reader.lines() {
                        if line.is_ok() {
                            total_errors += 1;
                        }
                    }
                }
            }
        }
    }

    info!("Total base errors: {}", total_errors);
    Ok(total_errors)
}
