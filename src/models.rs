/// Common data structures and models used across the codebase

use std::collections::HashMap;

/// Represents structural variant info
#[derive(Debug, Clone)]
pub struct StructuralVariant {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub sv_type: SVType,
    pub size: i64,
    pub support: usize,
    pub genotype: Option<Genotype>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SVType {
    Insertion,
    Deletion,
    Inversion,
    Translocation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    Homozygous,
    Heterozygous,
    Unknown,
}

/// Represents a base-level error (SNP or small indel)
#[derive(Debug, Clone)]
pub struct BaseError {
    pub chrom: String,
    pub position: u64,
    pub ref_base: char,
    pub alt_base: String,
    pub error_type: BaseErrorType,
    pub support: usize,
    pub depth: usize,
    pub p_value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseErrorType {
    SNP,
    Insertion,
    Deletion,
}

/// Contig metadata
#[derive(Debug, Clone)]
pub struct ContigInfo {
    pub name: String,
    pub length: u64,
    pub gc_content: Option<f64>,
    pub n50_rank: Option<usize>,
}

/// Alignment record (subset of SAM/BAM fields relevant to Inspector)
#[derive(Debug, Clone)]
pub struct AlignmentRecord {
    pub query_name: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u64,
    pub mapq: u8,
    pub cigar: String,
    pub query_len: u32,
    pub seq: Vec<u8>,
}

/// Coverage statistics at a genomic position
#[derive(Debug, Clone)]
pub struct CoverageStats {
    pub position: u64,
    pub depth: u32,
    pub bases: HashMap<char, u32>, // A, C, G, T, N counts
}

/// Pileup information from samtools mpileup
#[derive(Debug, Clone)]
pub struct PileupInfo {
    pub chrom: String,
    pub position: u64,
    pub ref_base: char,
    pub depth: u32,
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
}
