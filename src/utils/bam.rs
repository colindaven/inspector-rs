/// BAM/SAM file utilities via subprocess
/// Uses samtools commands for compatibility without system dependencies

use std::process::{Command, Stdio};
use std::io::BufRead;
use anyhow::{Result, Context};

/// Parse SAM/BAM file line (simplified - key fields only)
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u64,
    pub mapq: u8,
    pub cigar: String,
    pub seq_len: usize,
    pub nm_tag: Option<u32>, // edit distance
}

impl SamRecord {
    /// Parse a SAM line (tab-delimited)
    pub fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return None;
        }

        let flag = fields[1].parse().ok()?;
        let mapq = fields[4].parse().ok()?;
        let pos = fields[3].parse().ok()?;

        // Extract NM tag from optional fields
        let mut nm_tag = None;
        for i in 11..fields.len() {
            if fields[i].starts_with("NM:i:") {
                nm_tag = fields[i][5..].parse().ok();
                break;
            }
        }

        Some(SamRecord {
            qname: fields[0].to_string(),
            flag,
            rname: fields[2].to_string(),
            pos,
            mapq,
            cigar: fields[5].to_string(),
            seq_len: fields[9].len(),
            nm_tag,
        })
    }
}

/// Get effective mapping quality threshold
pub fn get_mapq_threshold(datatype: &str) -> u8 {
    match datatype {
        "hifi" => 10,
        "nanopore_94" | "nanopore_1041" => 5,
        _ => 0, // clr or others
    }
}

/// Get default minimum identity threshold
pub fn get_min_identity(datatype: &str) -> f64 {
    match datatype {
        "hifi" => 0.99,
        "nanopore_94" | "nanopore_1041" => 0.85,
        _ => 0.80, // clr or others
    }
}

/// Run samtools command and get output
pub fn run_samtools(args: &[&str]) -> Result<String> {
    let output = Command::new("samtools")
        .args(args)
        .stdout(Stdio::piped())
        .output()
        .context("Failed to execute samtools")?;

    if !output.status.success() {
        anyhow::bail!("samtools failed: {}", String::from_utf8_lossy(&output.stderr));
    }

    Ok(String::from_utf8(output.stdout)?)
}

/// Get BAM file statistics (number of mapped/unmapped reads)
pub fn get_bam_stats(bam_file: &str) -> Result<(u64, u64)> {
    let output = run_samtools(&["flagstat", bam_file])?;
    
    let mut total = 0u64;
    let mut mapped = 0u64;

    for line in output.lines() {
        if line.contains("mapped") && !line.contains("primary") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if let Ok(count) = parts[0].parse::<u64>() {
                mapped = count;
            }
        }
        if line.contains("in total") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if let Ok(count) = parts[0].parse::<u64>() {
                total = count;
            }
        }
    }

    Ok((total, mapped))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sam_record_parse() {
        let line = "read1\t0\tchr1\t1000\t20\t100M\t*\t0\t0\tACGT\t****\tNM:i:2";
        let record = SamRecord::from_line(line).unwrap();
        assert_eq!(record.qname, "read1");
        assert_eq!(record.flag, 0);
        assert_eq!(record.rname, "chr1");
        assert_eq!(record.pos, 1000);
        assert_eq!(record.mapq, 20);
        assert_eq!(record.cigar, "100M");
        assert_eq!(record.nm_tag, Some(2));
    }

    #[test]
    fn test_mapq_threshold() {
        assert_eq!(get_mapq_threshold("hifi"), 10);
        assert_eq!(get_mapq_threshold("nanopore_94"), 5);
        assert_eq!(get_mapq_threshold("nanopore_1041"), 5);
        assert_eq!(get_mapq_threshold("clr"), 0);
    }
}
