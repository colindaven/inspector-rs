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
                // Do NOT add to query_pos: query_length should equal the aligned
                // portion only (M+I+X), matching Python's query_alignment_length.
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

/// Parse CIGAR from split-read alignment.
/// Matches Python's `cigardeletion()`: captures indels >= 5bp, merges nearby
/// ones within the read, then retains only those >= min_size.
pub fn parse_cigar(
    _flag: u16,
    _chrom: &str,
    position: u64,
    cigar: &str,
    min_size: usize,
    max_size: usize,
) -> CigarInfo {
    // Step 1: parse with capture_min=5 to get small indels that may merge
    let capture_min = 5usize;
    let mut info = parse_cigar_ref(_flag, _chrom, position, cigar, capture_min, max_size);
    // Step 2: merge nearby indels within the read (Python's cigardeletion merge logic)
    merge_indels(&mut info.insertions, &mut info.deletions);
    // Step 3: filter to actual min_size
    info.insertions.retain(|&(_, len)| len >= min_size as u64);
    info.deletions.retain(|&(_, len)| len >= min_size as u64);
    info
}

/// Merge nearby indels within a read, mirroring Python's `cigardeletion()`.
///
/// Deletion merge: fixed 500bp gap window.
/// Insertion merge: adaptive window based on sizes of the pair being considered:
///   - max(l1,l2) < 100  → 600bp
///   - 100 ≤ max(l1,l2) < 500 → 400bp
///   - max(l1,l2) ≥ 500  → 600bp
///
/// Iterates from the rightmost pair backwards until no more merges are possible
/// (matching Python's single-pass-from-right with restart behaviour).
pub fn merge_indels(insertions: &mut Vec<(u64, u64)>, deletions: &mut Vec<(u64, u64)>) {
    // Merge deletions within fixed 500bp gap window
    let mut merged = true;
    while merged {
        merged = false;
        if deletions.len() <= 1 { break; }
        deletions.sort_by_key(|x| x.0);
        let mut i = deletions.len() - 1;
        while i > 0 {
            let gap = deletions[i].0.saturating_sub(deletions[i - 1].0 + deletions[i - 1].1);
            if gap <= 500 {
                let combined = (deletions[i - 1].0, deletions[i - 1].1 + deletions[i].1);
                deletions.remove(i);
                deletions[i - 1] = combined;
                merged = true;
                break;
            }
            i -= 1;
        }
    }

    // Merge insertions with adaptive window (matches Python exactly)
    merged = true;
    while merged {
        merged = false;
        if insertions.len() <= 1 { break; }
        insertions.sort_by_key(|x| x.0);
        let mut i = insertions.len() - 1;
        while i > 0 {
            let l1 = insertions[i].1;
            let l2 = insertions[i - 1].1;
            let max_len = l1.max(l2);
            // Python logic:
            //   window=200 if max(l1,l2)<100 else 400
            //   window=400 if window==400 and max(l1,l2)<500 else 600
            // Result: <100 → 600, 100..499 → 400, >=500 → 600
            let window: u64 = if max_len < 100 {
                600
            } else if max_len < 500 {
                400
            } else {
                600
            };
            let gap = insertions[i].0.saturating_sub(insertions[i - 1].0);
            if gap <= window {
                let combined = (insertions[i - 1].0, l1 + l2);
                insertions.remove(i);
                insertions[i - 1] = combined;
                merged = true;
                break;
            }
            i -= 1;
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
        let info = parse_cigar(0, "chr1", 1000, cigar, 5, 4000000);
        assert_eq!(info.left_clip, 10);
        assert_eq!(info.right_clip, 10);
        assert_eq!(info.query_length, 50); // aligned portion only (M), clips excluded
    }
}

