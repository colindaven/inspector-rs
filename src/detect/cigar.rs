/// CIGAR string parsing for identifying structural variations
/// Efficient byte-by-byte parsing without regex

use crate::models::StructuralVariant;

/// A single CIGAR operation (count, operation type)
#[derive(Debug, Clone, Copy)]
struct CigarOp {
    count: u64,
    op: char,
}

/// Parse raw CIGAR string into operations
/// Example: "100M2I100M" -> [(100, 'M'), (2, 'I'), (100, 'M')]
fn parse_cigar_ops(cigar: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num_str = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_str.push(c);
        } else {
            if let Ok(count) = num_str.parse::<u64>() {
                ops.push(CigarOp { count, op: c });
            }
            num_str.clear();
        }
    }
    ops
}

/// Result of CIGAR parsing
#[derive(Debug, Clone)]
pub struct CigarInfo {
    pub insertions: Vec<(u64, u64)>, // (position, length)
    pub deletions: Vec<(u64, u64)>,  // (position, length)
    pub ref_length: u64,
    pub query_length: u64,
    pub left_clip: u64,
    pub right_clip: u64,
}

/// Parse CIGAR from reference-to-contig alignment
/// Used for identifying errors when comparing assembly to long reads
pub fn parse_cigar_ref(
    _flag: u16,
    _chrom: &str,
    position: u64,
    cigar: &str,
    min_size: usize,
    max_size: usize,
) -> CigarInfo {
    let ops = parse_cigar_ops(cigar);
    let mut ref_pos = position;
    let mut query_pos = 0u64;
    let mut insertions = Vec::new();
    let mut deletions = Vec::new();
    let mut left_clip = 0u64;
    let mut right_clip = 0u64;

    for (idx, op) in ops.iter().enumerate() {
        match op.op {
            'M' | '=' | 'X' => {
                ref_pos += op.count;
                query_pos += op.count;
            }
            'I' => {
                if op.count >= min_size as u64 && op.count <= max_size as u64 {
                    insertions.push((ref_pos, op.count));
                }
                query_pos += op.count;
            }
            'D' => {
                if op.count >= min_size as u64 && op.count <= max_size as u64 {
                    deletions.push((ref_pos, op.count));
                }
                ref_pos += op.count;
            }
            'N' | 'P' => {
                ref_pos += op.count;
            }
            'S' => {
                if idx == 0 {
                    left_clip = op.count;
                } else {
                    right_clip = op.count;
                }
                query_pos += op.count;
            }
            'H' => {
                // Hard clip - doesn't consume query
                if idx == 0 {
                    left_clip = op.count;
                } else {
                    right_clip = op.count;
                }
            }
            _ => {} // Unknown operations
        }
    }

    CigarInfo {
        insertions,
        deletions,
        ref_length: ref_pos - position,
        query_length: query_pos,
        left_clip,
        right_clip,
    }
}

/// Parse CIGAR from split-read alignment
/// Similar to ref parsing but tracks both directions
pub fn parse_cigar(
    _flag: u16,
    _chrom: &str,
    position: u64,
    cigar: &str,
    min_size: usize,
    max_size: usize,
) -> CigarInfo {
    parse_cigar_ref(_flag, _chrom, position, cigar, min_size, max_size)
}

/// Merge nearby indels within a window
/// Used to combine signals from overlapping CIGAR indels
pub fn merge_indels(insertions: &mut Vec<(u64, u64)>, deletions: &mut Vec<(u64, u64)>) {
    // Merge insertions within window
    let insertion_window = 200;
    let mut merged = true;
    while merged {
        merged = false;
        insertions.sort_by_key(|x| x.0);

        for i in (1..insertions.len()).rev() {
            if insertions[i].0 - insertions[i - 1].0 <= insertion_window {
                let combined = (insertions[i - 1].0, insertions[i - 1].1 + insertions[i].1);
                insertions.remove(i);
                insertions.remove(i - 1);
                insertions.push(combined);
                merged = true;
                break;
            }
        }
    }

    // Merge deletions within window
    let deletion_window = 500;
    merged = true;
    while merged {
        merged = false;
        deletions.sort_by_key(|x| x.0);

        for i in (1..deletions.len()).rev() {
            if deletions[i].0 - deletions[i - 1].0 <= deletion_window {
                let combined = (deletions[i - 1].0, deletions[i - 1].1 + deletions[i].1);
                deletions.remove(i);
                deletions.remove(i - 1);
                deletions.push(combined);
                merged = true;
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar_simple() {
        let cigar = "100M";
        let info = parse_cigar(0, "chr1", 1000, cigar, 50, 4000000);
        assert_eq!(info.ref_length, 100);
        assert_eq!(info.query_length, 100);
        assert_eq!(info.insertions.len(), 0);
        assert_eq!(info.deletions.len(), 0);
    }

    #[test]
    fn test_parse_cigar_with_indels() {
        let cigar = "50M100I50M";
        let info = parse_cigar(0, "chr1", 1000, cigar, 50, 4000000);
        assert_eq!(info.ref_length, 100); // 50 + 0 + 50
        assert_eq!(info.query_length, 200); // 50 + 100 + 50
        assert_eq!(info.insertions.len(), 1);
        assert_eq!(info.insertions[0].1, 100); // insertion length
    }

    #[test]
    fn test_parse_cigar_with_deletions() {
        let cigar = "50M100D50M";
        let info = parse_cigar(0, "chr1", 1000, cigar, 50, 4000000);
        assert_eq!(info.ref_length, 200); // 50 + 100 + 50
        assert_eq!(info.query_length, 100); // 50 + 0 + 50
        assert_eq!(info.deletions.len(), 1);
        assert_eq!(info.deletions[0].1, 100); // deletion length
    }

    #[test]
    fn test_parse_cigar_with_clipping() {
        let cigar = "10S50M10S";
        let info = parse_cigar(0, "chr1", 1000, cigar, 50, 4000000);
        assert_eq!(info.left_clip, 10);
        assert_eq!(info.right_clip, 10);
        assert_eq!(info.query_length, 70); // 10 + 50 + 10
    }
}

