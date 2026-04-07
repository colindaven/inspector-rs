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
    let mut all_lengths = Vec::new();
    let mut total_length = 0u64;
    let mut total_length_large = 0u64;

    // Write valid contigs and compute statistics
    let valid_fa_path = format!("{}valid_contig.fa", outpath);
    let mut valid_fa = File::create(&valid_fa_path)?;
    let contig_info_path = format!("{}contig_length_info", outpath);
    let mut contig_info = File::create(&contig_info_path)?;

    for contig in all_contigs {
        let seq_len = contig.sequence.len() as u64;
        all_lengths.push(seq_len);
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
            }

            // Write valid contig to FASTA
            writeln!(valid_fa, ">{}", contig.name)?;
            writeln!(valid_fa, "{}", contig.sequence)?;
        }
    }

    // Compute N50
    let n50 = utils::compute_n50(all_lengths);

    // Get largest contig length
    let largest_contig_length = contig_lengths.values().max().copied().unwrap_or(0);

    let stats = FastaStats {
        chromosomes: chromosomes.clone(),
        chromosomes_map,
        chromosomes_large: chromosomes_large.clone(),
        total_length,
        total_length_large,
        n50,
        largest_contig_length,
        contig_lengths,
    };

    // Write summary statistics
    let summary_path = format!("{}summary_statistics", outpath);
    let mut summary = File::create(&summary_path)?;
    writeln!(summary, "Contig number\t{}", stats.chromosomes.len())?;
    writeln!(summary, "Total contig length\t{}", stats.total_length)?;
    writeln!(summary, "Longest contig\t{}", stats.largest_contig_length)?;
    writeln!(summary, "N50\t{}", stats.n50)?;
    writeln!(summary, "N50 (large contigs)\t{}", 
        utils::compute_n50(
            stats.chromosomes_large.iter()
                .filter_map(|c| stats.contig_lengths.get(c).copied())
                .collect()
        )
    )?;

    info!("Contigs: {}", stats.chromosomes.len());
    info!("Total length: {}", utils::format_size(stats.total_length));
    info!("N50: {}", utils::format_size(stats.n50));
    info!("Large contigs (≥{}bp): {}", min_size_assembly_error, stats.chromosomes_large.len());

    Ok(stats)
}

/// Aggregate mapping statistics from depth files
pub fn mapping_info_contig(
    outpath: &str,
    large_chroms: &[String],
    small_chroms: &[String],
    total_length: u64,
    total_length_large: u64,
) -> Result<usize> {
    // Returns average coverage
    // TODO: Implement coverage calculation from depth files
    Ok(20)
}

/// Merge and report assembly errors from clustering
pub fn assembly_info_cluster(
    outpath: &str,
    min_size: usize,
    max_size: usize,
) -> Result<()> {
    // Placeholder
    Ok(())
}

/// Reference-based assembly error reporting
pub fn assembly_info_ref(outpath: &str) -> Result<()> {
    // Placeholder
    Ok(())
}

/// Get reference alignment info and coverage
pub fn get_ref_align_info(
    path: &str,
    total_length: u64,
) -> Result<()> {
    // Placeholder
    Ok(())
}

/// Count base-pair errors from reference alignment
pub fn basepair_error_ref(
    outpath: &str,
    largest_chr: &str,
) -> Result<()> {
    // Placeholder
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
    pub contig_lengths: HashMap<String, u64>,
}
