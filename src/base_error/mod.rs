/// Base-level error detection (SNPs and small indels)
/// Parses pileup and detects errors using statistical tests

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};

use anyhow::{Context, Result};
use log::{debug, info};
use regex::Regex;
use statrs::distribution::{Binomial, DiscreteCDF};

#[derive(Debug, Clone)]
struct BaseErrorCandidate {
    chrom: String,
    start: usize,
    end: usize,
    error_type: String,
    size: usize,
    support: usize,
    depth: usize,
}

fn parse_indel_sequences(info: &str, marker: char) -> Vec<String> {
    let mut sequences = Vec::new();
    let bytes = info.as_bytes();
    let mut index = 0;

    while index < bytes.len() {
        if bytes[index] != marker as u8 {
            index += 1;
            continue;
        }

        index += 1;
        let length_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }

        if index == length_start {
            continue;
        }

        let Ok(length) = info[length_start..index].parse::<usize>() else {
            continue;
        };

        let seq_end = index.saturating_add(length).min(bytes.len());
        sequences.push(info[index..seq_end].to_string());
        index = seq_end;
    }

    sequences
}

fn most_frequent_sequence(sequences: &[String]) -> Option<&str> {
    let mut best_sequence = None;
    let mut best_count = 0usize;

    for sequence in sequences {
        let count = sequences.iter().filter(|candidate| *candidate == sequence).count();
        if count > best_count {
            best_count = count;
            best_sequence = Some(sequence.as_str());
        }
    }

    best_sequence
}

fn normalize_pileup_info(bases_str: &str, start_re: &Regex) -> String {
    start_re
        .replace_all(bases_str, "")
        .chars()
        .map(|c| match c {
            ',' => '.',
            'a' | 'A' => 'A',
            't' | 'T' => 'T',
            'c' | 'C' => 'C',
            'g' | 'G' => 'G',
            other => other,
        })
        .collect()
}

fn write_candidate(output: &mut File, candidate: &BaseErrorCandidate) -> Result<()> {
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        candidate.chrom,
        candidate.start,
        candidate.end,
        candidate.error_type,
        candidate.size,
        candidate.support,
        candidate.depth,
    )?;
    Ok(())
}

fn candidate_p_value(support: usize, depth: usize, propvalue: f64) -> Result<f64> {
    let distribution = Binomial::new(propvalue, depth as u64)
        .context("failed to construct binomial distribution")?;
    let threshold = support.saturating_sub(1) as u64;
    Ok(1.0 - distribution.cdf(threshold))
}

/// Detect SNVs and small indels from pileup
pub fn get_snv(
    path: &str,
    chrom: &str,
    contig_length: u64,
    mincount: usize,
    maxcov: usize,
    mindepth: Option<usize>,
    min_mapping_quality: u8,
) -> Result<()> {
    //info!("Detecting SNVs/indels for {}", chrom);

    let bam_file = format!("{}read_to_contig.bam", path);
    let fa_file = format!("{}valid_contig.fa", path);
    let workspace = format!("{}base_error_workspace/", path);
    let output_file = format!("{}baseerror_{}.bed", workspace, chrom);

    let indent_re = Regex::new(r"\+\d+[ACGTNacgtn]+")
        .context("indent_re compile failed")?;
    let del_re = Regex::new(r"-\d+[ACGTNacgtn]+")
        .context("del_re compile failed")?;
    let start_re = Regex::new(r"\^.")
        .context("start_re compile failed")?;

    let mut mpileup = Command::new("samtools")
        .arg("mpileup")
        .arg("-Q").arg("0")
        .arg("-q").arg(min_mapping_quality.to_string())
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
        let Ok(raw_depth) = fields[3].parse::<usize>() else {
            continue;
        };

        if ref_base != "N" && raw_depth >= mindepth_val && raw_depth <= maxcov {
            valid_bases += 1;
        }

        let info = normalize_pileup_info(fields[4], &start_re);
        let deleted_bases = info.matches('*').count();
        let effective_depth = raw_depth.saturating_sub(deleted_bases);

        if raw_depth < mincount || effective_depth > maxcov || effective_depth == 0 {
            continue;
        }

        let min_support = (effective_depth as f64 * 0.2).max(mincount as f64);
        let pos: usize = fields[1].parse().unwrap_or(0);
        let bed_start = pos.saturating_sub(1);
        let max_indel_len = mincount / 2;

        let insert_sequences = parse_indel_sequences(&info, '+');
        let delete_sequences = parse_indel_sequences(&info, '-');
        let eligible_insert_sequences: Vec<String> = insert_sequences
            .into_iter()
            .filter(|sequence| sequence.len() <= max_indel_len)
            .collect();
        let eligible_delete_sequences: Vec<String> = delete_sequences
            .into_iter()
            .filter(|sequence| sequence.len() <= max_indel_len)
            .collect();

        if (eligible_insert_sequences.len() as f64) >= min_support {
            if let Some(sequence) = most_frequent_sequence(&eligible_insert_sequences) {
                debug!("Found insertion at {}", fields[1]);
                error_count += 1;
                write_candidate(&mut output, &BaseErrorCandidate {
                    chrom: fields[0].to_string(),
                    start: bed_start,
                    end: pos,
                    error_type: "insertion".to_string(),
                    size: sequence.len(),
                    support: eligible_insert_sequences.len(),
                    depth: effective_depth,
                })?;
            }
        }

        if (eligible_delete_sequences.len() as f64) >= min_support {
            if let Some(sequence) = most_frequent_sequence(&eligible_delete_sequences) {
                debug!("Found deletion at {}", fields[1]);
                error_count += 1;
                let deletion_end = pos.saturating_add(sequence.len()).saturating_sub(1);
                write_candidate(&mut output, &BaseErrorCandidate {
                    chrom: fields[0].to_string(),
                    start: bed_start,
                    end: deletion_end,
                    error_type: "deletion".to_string(),
                    size: sequence.len(),
                    support: eligible_delete_sequences.len(),
                    depth: effective_depth,
                })?;
            }
        }

        if (info.matches('.').count() + deleted_bases) as f64 > 0.8 * raw_depth as f64 {
            continue;
        }

        let bases_no_insertions = indent_re.replace_all(&info, "");
        let bases = del_re.replace_all(&bases_no_insertions, "");
        let acount = bases.matches('A').count();
        let tcount = bases.matches('T').count();
        let ccount = bases.matches('C').count();
        let gcount = bases.matches('G').count();
        let candidates = [
            ('A', acount),
            ('T', tcount),
            ('C', ccount),
            ('G', gcount),
        ];

        if let Some((_, count)) = candidates.iter().max_by_key(|(_, count)| *count) {
            if (*count as f64) >= min_support {
                debug!("Found SNP at {}", fields[1]);
                error_count += 1;
                write_candidate(&mut output, &BaseErrorCandidate {
                    chrom: fields[0].to_string(),
                    start: bed_start,
                    end: pos,
                    error_type: "SNP".to_string(),
                    size: 1,
                    support: *count,
                    depth: effective_depth,
                })?;
            }
        }
    }

    let status = mpileup.wait()?;
    if !status.success() {
        anyhow::bail!(
            "samtools mpileup failed for contig '{}'\n  BAM: {}\n  Reference: {}",
            chrom,
            bam_file,
            fa_file
        );
    }

    let validbase_file = format!("{}validbase_{}", workspace, chrom);
    let mut validbase = File::create(&validbase_file)?;
    writeln!(validbase, "{}", valid_bases)?;

    let raw_error_percent = if contig_length == 0 {
        0.0
    } else {
        (error_count as f64 / contig_length as f64) * 100.0
    };
    // maxcov = contig_coverage * 2, so contig_coverage = maxcov / 2
    let contig_coverage = maxcov / 2;
    info!(
        "Found {} raw base-error candidates in {} ({:.4}% of {} bp, coverage={}x)",
        error_count,
        chrom,
        raw_error_percent,
        contig_length,
        contig_coverage
    );
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
    ctg_total_length: u64,
    datatype: &str,
    ave_depth: usize,
) -> Result<BaseErrorCounts> {
    debug!("Counting base errors");

    let workspace = format!("{}base_error_workspace/", path);
    let mut counts = BaseErrorCounts::default();
    let filtered_file = format!("{}small_scale_error.bed", path);
    let mut filtered_candidates = Vec::new();

    let (propvalue, pcutoff, readcutoff) = if datatype == "hifi" {
        let pcutoff = if ave_depth < 15 {
            0.1
        } else if ave_depth < 25 {
            0.02
        } else {
            0.01
        };
        (0.5, pcutoff, 0.75)
    } else {
        let pcutoff = if ave_depth < 25 { 0.1 } else { 0.05 };
        (0.4, pcutoff, 0.5)
    };

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

    // Apply Python-style filtering to baseerror_*.bed files and write filtered output.
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
                        let fields: Vec<&str> = line.split('\t').collect();
                        if fields.len() < 7 {
                            continue;
                        }

                        let support = fields[5].parse::<usize>().unwrap_or(0);
                        let depth = fields[6].parse::<usize>().unwrap_or(0);
                        if depth == 0 || (support as f64) < readcutoff * depth as f64 {
                            continue;
                        }

                        let pvalue = candidate_p_value(support, depth, propvalue)?;
                        if pvalue >= pcutoff {
                            continue;
                        }

                        counts.total += 1;
                        let typ = fields[3];
                        match typ {
                            "insertion"  => counts.small_scale_collapse   += 1,
                            "deletion"   => counts.small_scale_expansion  += 1,
                            "SNP"        => counts.base_substitution       += 1,
                            _ => {}
                        }

                        filtered_candidates.push(format!("{}\t{:.12e}", line, pvalue));
                    }
                }
            }
        }
    }

    let mut filtered_output = File::create(&filtered_file)?;
    writeln!(
        filtered_output,
        "#Contig_Name\tStart_Position\tEnd_Position\tType\tSize\tSupporting_Read\tDepth\tPvalue"
    )?;
    for candidate in &filtered_candidates {
        writeln!(filtered_output, "{}", candidate)?;
    }

    if ctg_total_length > 0 {
        let per_mbp = counts.total as f64 / ctg_total_length as f64 * 1_000_000.0;
        info!("Small-scale assembly error /per Mbp: {}", per_mbp);
    }

    let filtered_error_percent = if counts.valid_bases == 0 {
        0.0
    } else {
        (counts.total as f64 / counts.valid_bases as f64) * 100.0
    };

    info!(
        "Base errors — total: {} ({:.4}% of {} valid bases), substitution: {}, expansion: {}, collapse: {}",
        counts.total,
        filtered_error_percent,
        counts.valid_bases,
        counts.base_substitution,
        counts.small_scale_expansion,
        counts.small_scale_collapse,
    );
    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::{most_frequent_sequence, normalize_pileup_info, parse_indel_sequences};
    use regex::Regex;

    #[test]
    fn normalize_pileup_info_matches_python_cleanup() {
        let start_re = Regex::new(r"\^.").unwrap();
        let normalized = normalize_pileup_info(".,aA^]tTcCgG", &start_re);
        assert_eq!(normalized, "..AATTCCGG");
    }

    #[test]
    fn parse_indel_sequences_extracts_each_event() {
        let info = ".+5AACGT-2TT+1A";
        assert_eq!(parse_indel_sequences(info, '+'), vec!["AACGT", "A"]);
        assert_eq!(parse_indel_sequences(info, '-'), vec!["TT"]);
    }

    #[test]
    fn long_indels_are_filtered_before_support_count() {
        let parsed = parse_indel_sequences(".+5AACGT+2TT+1A", '+');
        let eligible: Vec<String> = parsed
            .into_iter()
            .filter(|sequence| sequence.len() <= 2)
            .collect();

        assert_eq!(eligible, vec!["TT", "A"]);
        assert_eq!(most_frequent_sequence(&eligible), Some("TT"));
    }
}
