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

    // Compile regexes ONCE — not per-line — for speed
    let indent_re = Regex::new(r"\+\d+[ACGTNacgtn]+").context("indent_re compile failed")?;
    let del_re    = Regex::new(r"-\d+[ACGTNacgtn]+").context("del_re compile failed")?;
    let start_re  = Regex::new(r"\^.").context("start_re compile failed")?;

    // Stream mpileup DIRECTLY via piped stdout — avoids writing a .pileup file to disk
    // and re-reading it, cutting I/O in half per contig.
    let mut mpileup = Command::new("samtools")
        .arg("mpileup")
        .arg("-Q").arg("0")
        .arg("-r").arg(chrom)
        .arg("-f").arg(&fa_file)
        .arg(&bam_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| anyhow::anyhow!(
            "Failed to run samtools mpileup: {}\n  BAM: {}\n  Reference: {}\n  Contig: {}\n  Hint: ensure samtools is in PATH and both files exist",
            e, bam_file, fa_file, chrom
        ))?;

    let reader = BufReader::new(mpileup.stdout.take().unwrap());
    let mut output = File::create(&output_file)?;
    use std::io::Write;

    let mindepth_val = mindepth.unwrap_or(maxcov / 10);
    let mut valid_bases = 0u64;
    let mut error_count = 0u64;

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

        // Normalise the base string in-place using pre-compiled regexes
        // Remove insertion/deletion sequence bodies before counting * chars
        let bases_no_ins  = indent_re.replace_all(bases_str, "");
        let bases_no_indels = del_re.replace_all(&*bases_no_ins, "");
        // Remove read-start markers (^<mapq-char>)
        let bases_no_start = start_re.replace_all(&*bases_no_indels, "");
        // Normalise case: lowercase = reverse-strand match, uppercase mismatch
        let bases: String = bases_no_start.chars()
            .map(|c| match c { ',' => '.', 'a'|'A' => 'A', 't'|'T' => 'T',
                               'c'|'C' => 'C', 'g'|'G' => 'G', other => other })
            .collect();

        let deletions = bases.matches('*').count();
        let effective_depth = depth.saturating_sub(deletions);

        if effective_depth == 0 {
            continue;
        }

        // Count raw insertion events (each + in the original string is one event)
        let insertions = indent_re.find_iter(bases_str).count();

        let min_support = ((depth as f64 * 0.2).max(mincount as f64)) as usize;

        let pos: usize = fields[1].parse().unwrap_or(0);
        let bed_start = pos.saturating_sub(1);

        if insertions >= min_support {
            debug!("Found insertion at {}", fields[1]);
            error_count += 1;
            writeln!(output, "{}\t{}\t{}\tinsertion\t{}", fields[0], bed_start, pos, insertions)?;
        }

        if deletions >= min_support {
            debug!("Found deletion at {}", fields[1]);
            error_count += 1;
            writeln!(output, "{}\t{}\t{}\tdeletion\t{}", fields[0], bed_start, pos, deletions)?;
        }

        // Detect SNPs: count non-reference, non-gap, non-marker bases
        let mut base_counts = std::collections::HashMap::new();
        for c in bases.chars() {
            if c != '.' && c != '-' && c != '+' && c != '*' && c != '^' {
                *base_counts.entry(c).or_insert(0usize) += 1;
            }
        }

        for (_alt_base, count) in base_counts {
            if count >= min_support && count < effective_depth {
                let freq = count as f64 / effective_depth as f64;
                if freq > 0.1 && freq < 0.9 {
                    debug!("Found SNP at {}", fields[1]);
                    error_count += 1;
                    writeln!(output, "{}\t{}\t{}\tSNP\t{}", fields[0], bed_start, pos, count)?;
                }
            }
        }
    }

    // Wait for mpileup to finish and propagate non-zero exit
    let status = mpileup.wait()?;
    if !status.success() {
        anyhow::bail!("samtools mpileup failed for contig '{}'\n  BAM: {}\n  Reference: {}",
            chrom, bam_file, fa_file);
    }

    // Write per-chrom validbase count (one file per contig — safe across parallel calls)
    let validbase_file = format!("{}validbase_{}", workspace, chrom);
    let mut validbase = File::create(&validbase_file)?;
    writeln!(validbase, "{}", valid_bases)?;

    info!("Found {} errors in {}", error_count, chrom);
    Ok(())
}

/// Breakdown of base-scale errors by type, matching Python's summary categories.
#[derive(Debug, Clone, Default)]
pub struct BaseErrorCounts {
    /// Total filtered errors
    pub total: usize,
    /// Insertions in pileup = Small-scale collapse (assembly missing bases)
    pub small_scale_collapse: usize,
    /// Deletions in pileup = Small-scale expansion (assembly has extra bases)
    pub small_scale_expansion: usize,
    /// SNP mismatches = Base substitution
    pub base_substitution: usize,
    /// Total valid (covered) bases used for QV denominator
    pub valid_bases: u64,
}

/// Count and categorise base errors from all per-chrom bed files.
pub fn count_base_errors(
    path: &str,
    _ctg_total_length: u64,
    _datatype: &str,
    _ave_depth: usize,
) -> Result<BaseErrorCounts> {
    debug!("Counting base errors");

    let workspace = format!("{}base_error_workspace/", path);
    let mut counts = BaseErrorCounts::default();

    // Sum valid_bases from all per-chrom validbase_{chrom} files
    if let Ok(entries) = std::fs::read_dir(&workspace) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.file_name().and_then(|n| n.to_str())
                .map(|n| n.starts_with("validbase_"))
                .unwrap_or(false)
            {
                if let Ok(s) = std::fs::read_to_string(&p) {
                    if let Ok(v) = s.trim().parse::<u64>() {
                        counts.valid_bases += v;
                    }
                }
            }
        }
    }

    // Count errors by type from baseerror_*.bed files
    if let Ok(entries) = std::fs::read_dir(&workspace) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.file_name().and_then(|n| n.to_str())
                .map(|n| n.starts_with("baseerror_") && n.ends_with(".bed"))
                .unwrap_or(false)
            {
                if let Ok(file) = File::open(&p) {
                    for line in BufReader::new(file).lines().flatten() {
                        if line.is_empty() { continue; }
                        counts.total += 1;
                        // Field 3 (0-indexed) is the type
                        let typ = line.split('\t').nth(3).unwrap_or("");
                        match typ {
                            "insertion"  => counts.small_scale_collapse   += 1,
                            "deletion"   => counts.small_scale_expansion  += 1,
                            "SNP"        => counts.base_substitution       += 1,
                            _ => {}
                        }
                    }
                }
            }
        }
    }

    info!("Base errors — total: {}, substitution: {}, expansion: {}, collapse: {}, valid_bases: {}",
        counts.total, counts.base_substitution, counts.small_scale_expansion,
        counts.small_scale_collapse, counts.valid_bases);
    Ok(counts)
}
