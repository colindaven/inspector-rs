/// SV detection from BAM alignments
/// Parses CIGAR strings, extracts split reads, and identifies structural error signals

pub mod cigar;
pub mod segment;

pub use cigar::{parse_cigar, parse_cigar_ref, merge_indels, CigarInfo};
pub use segment::{detect_segments, ReadSegment};

use std::process::{Command, Stdio};
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::collections::HashMap;
use anyhow::{Result, Context};
use log::{debug, info};
use crate::utils::bam::SamRecord;

/// Main entry point for SV detection on a single contig
pub fn detect_sortbam(
    workpath: &str,
    min_size: usize,
    max_size: usize,
    chrom: &str,
) -> Result<()> {
    info!("SV detection for contig: {}", chrom);

    let bam_file = format!("{}read_to_contig.bam", workpath);
    let output_file = format!("{}debreak_workspace/read_to_contig_{}.debreak.temp", workpath, chrom);

    // Use samtools view to stream BAM records as SAM
    let samtools = Command::new("samtools")
        .arg("view")
        .arg("-h")
        .arg(&bam_file)
        .arg(chrom)
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| anyhow::anyhow!("Failed to start samtools view: {}\n  BAM file: {}\n  Contig: {}\n  Hint: ensure samtools is in PATH and the BAM file is indexed", e, bam_file, chrom))?;

    let reader = BufReader::new(samtools.stdout.unwrap());
    let mut sv_calls = Vec::new();
    let mut segment_reads: HashMap<String, Vec<ReadSegment>> = HashMap::new();
    let mut total_map_length = 0u64;
    let mut number_reads = 0u64;

    for line in reader.lines() {
        let line = line?;
        
        // Skip header lines
        if line.starts_with('@') {
            continue;
        }

        let Some(record) = SamRecord::from_line(&line) else {
            continue;
        };

        // Skip secondary alignments (we only want primary + supplementary)
        if (record.flag & 256) != 0 {
            continue; // secondary
        }

        number_reads += 1;
        total_map_length += record.seq_len as u64;

        // Parse CIGAR to get indel signals
        let cigar_info = parse_cigar(record.flag, &record.rname, record.pos, &record.cigar, min_size, max_size);

        // Write CIGAR-based signals immediately
        for (pos, len) in &cigar_info.insertions {
            sv_calls.push(format!(
                "{}\t{}\t{}\tI-cigar\t{}",
                chrom, pos, len, record.qname
            ));
        }
        for (pos, len) in &cigar_info.deletions {
            sv_calls.push(format!(
                "{}\t{}\t{}\tD-cigar\t{}",
                chrom, pos, len, record.qname
            ));
        }

        // Collect split reads by name for segment analysis
        let segment = ReadSegment {
            query_name: record.qname.clone(),
            flag: record.flag,
            rname: record.rname.clone(),
            pos: record.pos,
            refend: record.pos + cigar_info.ref_length,
            cigar_info: (cigar_info.left_clip, cigar_info.query_length, cigar_info.right_clip),
            mapq: record.mapq,
        };

        segment_reads.entry(record.qname).or_insert_with(Vec::new).push(segment);
    }

    debug!("Processed {} reads with {} total mapped bp", number_reads, total_map_length);

    // Detect split-read based SVs
    for (_read_name, segments) in segment_reads {
        if segments.len() > 1 {
            let split_svs = detect_segments(&segments, min_size, max_size);
            sv_calls.extend(split_svs);
        }
    }

    // Write output
    let mut output = File::create(&output_file)?;
    use std::io::Write;
    for sv_call in sv_calls {
        writeln!(output, "{}", sv_call)?;
    }

    info!("Wrote {} SV calls to {}", number_reads, output_file);
    Ok(())
}

/// Lighter version for contigs where SV detection is skipped
pub fn detect_sortbam_nosv(
    workpath: &str,
    chrom: &str,
    _contig_type: &str,
) -> Result<()> {
    debug!("Collecting coverage for contig: {} (no SV detection)", chrom);

    let bam_file = format!("{}read_to_contig.bam", workpath);
    let output_file = format!("{}map_depth/read_to_contig_{}.depth", workpath, chrom);

    let samtools = Command::new("samtools")
        .arg("view")
        .arg(&bam_file)
        .arg(chrom)
        .stdout(Stdio::piped())
        .spawn()
        .context("Failed to start samtools view")?;

    let reader = BufReader::new(samtools.stdout.unwrap());
    let mut total_map_length = 0u64;
    let mut number_reads = 0u64;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let Some(record) = SamRecord::from_line(&line) else {
            continue;
        };

        if (record.flag & 256) == 0 { // Not secondary
            number_reads += 1;
            total_map_length += record.seq_len as u64;
        }
    }

    // Write depth summary
    let mut output = File::create(&output_file)?;
    use std::io::Write;
    writeln!(output, "{}\t{}\t{}", chrom, number_reads, total_map_length)?;

    Ok(())
}

/// Reference-based SV detection from SAM files
pub fn detect_sam_ref(
    _filename: &str,
    _readpath: &str,
    _writepath: &str,
    _min_size: usize,
    _max_size: usize,
) -> Result<()> {
    // Placeholder
    Ok(())
}
