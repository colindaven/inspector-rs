/// Segment-based SV detection from split reads
/// Groups supplementary alignments from same read to detect SVs

use std::collections::HashMap;
use log::debug;

/// Segment information from alignment
#[derive(Debug, Clone)]
pub struct ReadSegment {
    pub query_name: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u64,
    pub refend: u64,
    pub cigar_info: (u64, u64, u64), // (left_clip, aligned_len, right_clip)
    pub mapq: u8,
}

/// Detect structural variants from split-read segments
pub fn detect_segments(
    segments: &[ReadSegment],
    min_size: usize,
    max_size: usize,
) -> Vec<String> {
    let mut sv_calls = Vec::new();

    // Group segments by read name
    let mut read_groups: HashMap<String, Vec<ReadSegment>> = HashMap::new();
    for segment in segments {
        read_groups.entry(segment.query_name.clone())
            .or_insert_with(Vec::new)
            .push(segment.clone());
    }

    // Sort each read's segments by position
    for (_read_name, mut segs) in read_groups {
        if segs.len() < 2 {
            continue; // Need at least 2 segments for SV detection
        }

        segs.sort_by_key(|s| s.pos);

        // Check all pairs of segments for SVs
        for i in 0..segs.len() {
            for j in (i + 1)..segs.len() {
                let left = &segs[i];
                let right = &segs[j];

                // Same chromosome - insertions/deletions
                if left.rname == right.rname {
                    detect_insertion_deletion(left, right, min_size, max_size, &mut sv_calls);
                    detect_inversion(left, right, min_size, max_size, &mut sv_calls);
                } else {
                    // Different chromosomes - translocation
                    detect_translocation(left, right, &mut sv_calls);
                }
            }
        }
    }

    sv_calls
}

/// Detect insertions and deletions
fn detect_insertion_deletion(
    left: &ReadSegment,
    right: &ReadSegment,
    min_size: usize,
    max_size: usize,
    sv_calls: &mut Vec<String>,
) {
    let window_max = 500;
    let overlap_window = -200i64;

    // Calculate reference and query gaps between segments
    let ref_gap = right.pos as i64 - left.refend as i64;
    
    // Query gap: distance in query between end of left segment and start of right segment
    // If segments are consecutive on the query, query_gap ≈ 0
    // If there's a gap, query_gap > 0; if overlap, query_gap < 0
    let left_query_end = left.cigar_info.0 as i64 + left.cigar_info.1 as i64;
    let right_query_start = right.cigar_info.0 as i64;
    let query_gap = right_query_start - left_query_end;

    // Insertion: reference gap is small but query extends beyond
    // This means query bases align to the same ref position (insertion in assembly)
    if ref_gap > overlap_window && ref_gap < window_max as i64 && query_gap > 100 {
        let insert_size = query_gap;
        if insert_size >= min_size as i64 && insert_size <= max_size as i64 {
            let sv_line = format!(
                "{}\t{}\t{}\tINS-segment\t{}\t{}",
                left.rname, left.refend, insert_size,
                left.query_name, (left.mapq as u32 + right.mapq as u32) / 2
            );
            sv_calls.push(sv_line);
            return;
        }
    }

    // Deletion: reference gap is large compared to query gap
    // This means assembly is missing bases that are in the reads
    if ref_gap > 100 {
        let delete_size = ref_gap - query_gap;
        if delete_size >= min_size as i64 && delete_size <= max_size as i64 {
            let sv_line = format!(
                "{}\t{}\t{}\tDEL-segment\t{}\t{}",
                left.rname, left.refend, delete_size,
                left.query_name, (left.mapq as u32 + right.mapq as u32) / 2
            );
            sv_calls.push(sv_line);
        }
    }
}

/// Detect inversions
fn detect_inversion(
    left: &ReadSegment,
    right: &ReadSegment,
    min_size: usize,
    max_size: usize,
    sv_calls: &mut Vec<String>,
) {
    if left.rname != right.rname {
        return;
    }

    // Check orientation: left is forward, right is reverse = inversion
    let left_forward = (left.flag & 16) == 0;
    let right_forward = (right.flag & 16) == 0;

    if left_forward == right_forward {
        return; // Same orientation, not an inversion
    }

    let inv_size = (right.pos as i64 - left.refend as i64) as usize;
    if inv_size >= min_size && inv_size <= max_size {
        let sv_line = format!(
            "{}\t{}\t{}\tINV-segment\t{}\t{}",
            left.rname, left.refend, inv_size,
            left.query_name, (left.mapq as u32 + right.mapq as u32) / 2
        );
        sv_calls.push(sv_line);
    }
}

/// Detect translocations (different chromosomes)
fn detect_translocation(
    left: &ReadSegment,
    _right: &ReadSegment,
    sv_calls: &mut Vec<String>,
) {
    // Placeholder: translocations require more complex handling
    // involving inter-chromosomal alignment patterns
    let _sv_line = format!(
        "{}\t{}\tTRA\t{}\t{}",
        left.rname, left.pos, left.query_name, left.mapq
    );
    // sv_calls.push(sv_line); // TODO: Implement full translocation logic
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_segments_simple() {
        // Two segments on same read at different positions
        // This should detect a deletion (200bp gap on reference, 0bp on query)
        let segments = vec![
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 0,
                rname: "chr1".to_string(),
                pos: 1000,
                refend: 1100,
                cigar_info: (0, 100, 0), // (left_clip, aligned_len, right_clip)
                mapq: 20,
            },
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 0,
                rname: "chr1".to_string(),
                pos: 1300,
                refend: 1400,
                cigar_info: (0, 100, 0),
                mapq: 20,
            },
        ];

        let svs = detect_segments(&segments, 50, 4000000);
        // ref_gap = 1300-1100 = 200, query_gap = 0
        // delete_size = 200 - 0 = 200 >= 50 -> should detect
        assert!(!svs.is_empty(), "Expected deletion of 200bp, got {:?}", svs);
        assert!(svs[0].contains("DEL-segment"), "Expected DEL-segment call");
    }
}
