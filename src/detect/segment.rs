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
        if segs.len() < 2 || segs.len() > 20 {
            continue; // Need 2–20 segments; >20 is noise
        }

        // Require at least one forward-strand (flag bit 4 clear) alignment
        let has_forward = segs.iter().any(|s| (s.flag & 16) == 0);
        if !has_forward {
            continue;
        }

        // Drop segments with very short aligned length (< 300 bp) — mirrors Python's
        // supplementary-alignment quality gate in segmentdeletion()
        segs.retain(|s| s.cigar_info.1 >= 300);
        if segs.len() < 2 {
            continue;
        }

        segs.sort_by_key(|s| s.pos);

        // Check all pairs of segments for SVs
        for i in 0..segs.len() {
            for j in (i + 1)..segs.len() {
                let left = &segs[i];
                let right = &segs[j];

                // Same chromosome only
                if left.rname != right.rname {
                    continue;
                }

                let left_rev = (left.flag & 16) != 0;
                let right_rev = (right.flag & 16) != 0;

                if left_rev == right_rev {
                    // Same direction → insertions/deletions
                    // Python also requires: right.refend - left.refend > -200
                    if (right.refend as i64 - left.refend as i64) > -200 {
                        detect_insertion_deletion(left, right, min_size, max_size, &mut sv_calls);
                    }
                } else {
                    // Opposite direction → inversions
                    detect_inversion(left, right, min_size, max_size, &mut sv_calls);
                }
            }
        }
    }

    sv_calls
}

/// Detect insertions and deletions from same-direction split-read pairs.
/// Mirrors Python's `segmentdeletion` for `samedirchr` pairs.
///
/// Insertion formula (Python):
///   window = 300; if abs(right.pos - left.refend) <= window:
///     overlap = right.pos - left.refend
///     ins_size = right.leftclip - left.aligned - left.leftclip - overlap
///   Which simplifies to: ins_size = query_gap - ref_gap
///
/// Deletion formula (Python):
///   overlapmap = left.leftclip + left.aligned - right.leftclip  (= -query_gap)
///   if -200 < overlapmap < 2000:
///     del_size = (right.pos - left.refend) + overlapmap = ref_gap - query_gap
fn detect_insertion_deletion(
    left: &ReadSegment,
    right: &ReadSegment,
    min_size: usize,
    max_size: usize,
    sv_calls: &mut Vec<String>,
) {
    // ref_gap = right.pos - left.refend (Python's "overlap" or basis for it)
    let ref_gap = right.pos as i64 - left.refend as i64;

    // query_gap = right.leftclip - (left.leftclip + left.aligned)
    let left_query_end = left.cigar_info.0 as i64 + left.cigar_info.1 as i64;
    let right_query_start = right.cigar_info.0 as i64;
    let query_gap = right_query_start - left_query_end;

    // --- Insertion detection ---
    // Python: if abs(rightread[3]-leftread[4]) <= 300 (window=300)
    let ins_window: i64 = 300;
    if ref_gap.abs() <= ins_window {
        let insert_size = query_gap - ref_gap; // Python: rightinfo[0]-leftinfo[1]-leftinfo[0]-overlap
        if insert_size >= min_size as i64 && insert_size <= max_size as i64 {
            let pos = std::cmp::min(right.pos, left.refend);
            sv_calls.push(format!(
                "{}\t{}\t{}\tI-segment\t{}\t{}\t{}",
                left.rname, pos, insert_size,
                left.query_name,
                left.flag as u32 + right.flag as u32,
                (left.mapq as u32 + right.mapq as u32) / 2
            ));
            return;
        }
    }

    // --- Deletion detection ---
    // Python: overlapmap = leftinfo[0]+leftinfo[1]-rightinfo[0] = -query_gap
    //         if -200 < overlapmap < 2000:
    //           del_size = rightread[3]-leftread[4]+overlapmap = ref_gap + (-query_gap) = ref_gap - query_gap
    let overlapmap = -query_gap;
    let del_overlap_window: i64 = -200;
    let del_window_max: i64 = 2000;
    if overlapmap > del_overlap_window && overlapmap < del_window_max {
        let delete_size = ref_gap + overlapmap; // = ref_gap - query_gap
        if delete_size >= min_size as i64 && delete_size <= max_size as i64 {
            let pos = left.refend as i64 - std::cmp::max(0, overlapmap);
            sv_calls.push(format!(
                "{}\t{}\t{}\tD-segment\t{}\t{}\t{}",
                left.rname, pos.max(0), delete_size,
                left.query_name,
                left.flag as u32 + right.flag as u32,
                (left.mapq as u32 + right.mapq as u32) / 2
            ));
        }
    }
}

/// Detect inversions using Python's two-case overlapmap geometry check.
///
/// Python uses `flag % 32 > 15` to identify reverse-strand alignments (isolates
/// the lower 5 flag bits, checking bit 4 = 0x10 = reverse complement). We use
/// `flag & 16 != 0` for the same test.
///
/// Two orientation cases are tried (mirroring Python's segmentdeletion `samechr` loop):
///   Case 1 – inversion anchored at leftread's refend:
///     overlapmap = right.left_clip + right.aligned_len − left.right_clip
///     ref_span   = right.refend  − left.refend
///     pos        = left.refend
///   Case 2 – inversion anchored at leftread's pos:
///     overlapmap = right.aligned_len + right.right_clip − left.left_clip
///     ref_span   = right.pos − left.pos
///     pos        = left.pos
///
/// Acceptance gate (both cases): −200 < overlapmap < 500
///                                ref_span ≥ max(100, overlapmap)
///                                inv_size = ref_span − overlapmap ∈ [min_size, max_size]
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

    // Opposite strands required
    let left_rev = (left.flag & 16) != 0;
    let right_rev = (right.flag & 16) != 0;
    if left_rev == right_rev {
        return;
    }

    let overlap_window: i64 = -200;
    let window_max: i64 = 500;

    // cigar_info = (left_clip, aligned_len, right_clip)
    let (ll, la, lr) = (left.cigar_info.0 as i64, left.cigar_info.1 as i64, left.cigar_info.2 as i64);
    let (rl, ra, rr) = (right.cigar_info.0 as i64, right.cigar_info.1 as i64, right.cigar_info.2 as i64);

    // Case 1 — anchored at left.refend / right.refend
    {
        let overlapmap = rl + ra - lr;
        let ref_span   = right.refend as i64 - left.refend as i64;
        if overlapmap > overlap_window && overlapmap < window_max
            && ref_span >= std::cmp::max(100, overlapmap)
        {
            let inv_size = ref_span - overlapmap;
            if inv_size > 0 && inv_size as usize >= min_size && (inv_size as usize) <= max_size {
                debug!("INV-segment case1 {} pos={} size={}", left.rname, left.refend, inv_size);
                sv_calls.push(format!(
                    "{}\t{}\t{}\tINV-segment\t{}\t{}\t{}",
                    left.rname, left.refend, inv_size,
                    left.query_name,
                    left.flag as u32 + right.flag as u32,
                    (left.mapq as u32 + right.mapq as u32) / 2
                ));
                return; // one signal per pair is sufficient
            }
        }
    }

    // Case 2 — anchored at left.pos / right.pos
    {
        let overlapmap = ra + rr - ll;
        let ref_span   = right.pos as i64 - left.pos as i64;
        if overlapmap > overlap_window && overlapmap < window_max
            && ref_span >= std::cmp::max(100, overlapmap)
        {
            let inv_size = ref_span - overlapmap;
            if inv_size > 0 && inv_size as usize >= min_size && (inv_size as usize) <= max_size {
                debug!("INV-segment case2 {} pos={} size={}", left.rname, left.pos, inv_size);
                sv_calls.push(format!(
                    "{}\t{}\t{}\tINV-segment\t{}\t{}\t{}",
                    left.rname, left.pos, inv_size,
                    left.query_name,
                    left.flag as u32 + right.flag as u32,
                    (left.mapq as u32 + right.mapq as u32) / 2
                ));
            }
        }
    }
}

/// Detect translocations (different chromosomes)
fn detect_translocation(
    left: &ReadSegment,
    _right: &ReadSegment,
    _sv_calls: &mut Vec<String>,
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
        // cigar_info.1 must be >= 300 to pass the segment quality gate (mirrors Python's
        // 300bp minimum aligned-length check for supplementary alignments)
        let segments = vec![
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 0,
                rname: "chr1".to_string(),
                pos: 1000,
                refend: 1100,
                cigar_info: (0, 300, 0), // (left_clip, aligned_len, right_clip)
                mapq: 20,
            },
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 0,
                rname: "chr1".to_string(),
                pos: 1300,
                refend: 1600,
                cigar_info: (0, 300, 0),
                mapq: 20,
            },
        ];

        let svs = detect_segments(&segments, 50, 4000000);
        // ref_gap = 1300-1100 = 200, query_gap = 0-300 = -300
        // delete_size = 200 - (-300) = 500 >= 50 -> should detect
        assert!(!svs.is_empty(), "Expected deletion of 200bp, got {:?}", svs);
        assert!(svs[0].contains("D-segment"), "Expected D-segment call");
            // Verify 7-field format: chrom\tpos\tsize\ttype\treadname\tflag_sum\tavg_mapq
            let fields: Vec<&str> = svs[0].split('\t').collect();
            assert_eq!(fields.len(), 7, "Expected 7 tab-delimited fields in debreak.temp format, got {}: {:?}", fields.len(), svs[0]);
            assert_eq!(fields[5], "0",  "Field[5] should be flag sum (0+0=0)");
            assert_eq!(fields[6], "20", "Field[6] should be avg mapq (20+20)/2=20");
    }

    #[test]
    fn test_detect_inversion_geometric_check() {
        // A forward-strand segment followed by a reverse-strand segment.
        // With Python's overlapmap check:
        //   Case 1: overlapmap = right.left_clip + right.aligned_len - left.right_clip
        //                       = 50 + 300 - 50 = 300
        //           ref_span   = right.refend - left.refend = 5500 - 1500 = 4000
        //           Gate: -200 < 300 < 500 ✓; ref_span(4000) >= max(100,300) ✓
        //           inv_size = 4000 - 300 = 3700 → reported
        let segments = vec![
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 0,   // forward strand
                rname: "chr1".to_string(),
                pos: 1000,
                refend: 1500,
                cigar_info: (50, 400, 50), // left_clip=50, aligned_len=400, right_clip=50
                mapq: 30,
            },
            ReadSegment {
                query_name: "read1".to_string(),
                flag: 16,  // reverse strand
                rname: "chr1".to_string(),
                pos: 1200,
                refend: 5500,
                cigar_info: (50, 300, 100),
                mapq: 30,
            },
        ];

        let svs = detect_segments(&segments, 50, 4000000);
        assert!(!svs.is_empty(), "Expected inversion signal, got {:?}", svs);
        assert!(svs.iter().any(|s| s.contains("INV-segment")), "Expected INV-segment: {:?}", svs);
    }

    #[test]
    fn test_detect_inversion_rejects_same_strand() {
        // Both forward – must not produce INV
        let segments = vec![
            ReadSegment { query_name: "read1".to_string(), flag: 0, rname: "chr1".to_string(),
                pos: 1000, refend: 2000, cigar_info: (0, 1000, 0), mapq: 30 },
            ReadSegment { query_name: "read1".to_string(), flag: 0, rname: "chr1".to_string(),
                pos: 3000, refend: 4000, cigar_info: (0, 1000, 0), mapq: 30 },
        ];
        let svs = detect_segments(&segments, 50, 4000000);
        assert!(!svs.iter().any(|s| s.contains("INV-segment")),
            "Should not detect inversion for same-strand pair: {:?}", svs);
    }

    #[test]
    fn test_detect_inversion_rejects_short_segments() {
        // aligned_len = 100 (< 300) — must be filtered by the length pre-filter
        let segments = vec![
            ReadSegment { query_name: "read1".to_string(), flag: 0, rname: "chr1".to_string(),
                pos: 1000, refend: 1100, cigar_info: (0, 100, 0), mapq: 30 },
            ReadSegment { query_name: "read1".to_string(), flag: 16, rname: "chr1".to_string(),
                pos: 5000, refend: 5100, cigar_info: (0, 100, 0), mapq: 30 },
        ];
        let svs = detect_segments(&segments, 50, 4000000);
        assert!(!svs.iter().any(|s| s.contains("INV-segment")),
            "Should not detect inversion for short segments: {:?}", svs);
    }
}
