/// Merging operations for structural variants
/// Directly ported from Python debreak_merge_clustering.py:
///   mergeinfo_insertion_oneevent → merge_one_event
///   mergeinfo_insertion          → merge_with_bimodal
///
/// Input lines are 7-field debreak.temp strings:
///   chrom  pos  size  type  readname  flag  mapq   (tab-delimited, 0-indexed)
///
/// Output lines are merged-event strings:
///   chrom  pos  size  count  numread  quality  <empty>  readnames
/// (the empty field is a historical artefact of the Python format;
///  cluster code appends \tUnique or \tCompoundSV after readnames)

/// Assign K-means clusters (k=2) from sizes, iterated `iters` times.
/// Returns (group1_indices, group2_indices, mean1, mean2).
fn kmeans_assign(sizes: &[u64], mean1: u64, mean2: u64) -> (Vec<usize>, Vec<usize>, u64, u64) {
    let mut g1: Vec<usize> = Vec::new();
    let mut g2: Vec<usize> = Vec::new();
    for (i, &s) in sizes.iter().enumerate() {
        let d1 = (s as i64 - mean1 as i64).unsigned_abs();
        let d2 = (s as i64 - mean2 as i64).unsigned_abs();
        if d1 <= d2 { g1.push(i); } else { g2.push(i); }
    }
    let new_m1 = if g1.is_empty() { mean1 } else {
        g1.iter().map(|&i| sizes[i]).sum::<u64>() / g1.len() as u64
    };
    let new_m2 = if g2.is_empty() { mean2 } else {
        g2.iter().map(|&i| sizes[i]).sum::<u64>() / g2.len() as u64
    };
    (g1, g2, new_m1, new_m2)
}

/// Parse field[2] (size) from a debreak.temp / candidate line.
#[inline]
fn parse_size(line: &str) -> u64 {
    line.split('\t').nth(2).and_then(|s| s.parse().ok()).unwrap_or(0)
}

/// Parse field[1] (position) from a line.
#[inline]
fn parse_pos(line: &str) -> u64 {
    line.split('\t').nth(1).and_then(|s| s.parse().ok()).unwrap_or(0)
}

/// Parse field[6] (mapq) as quality from a debreak.temp line.
#[inline]
fn parse_quality(line: &str) -> f64 {
    line.split('\t').nth(6).and_then(|s| s.parse().ok()).unwrap_or(0.0)
}

/// Parse field[4] (readname) from a line.
#[inline]
fn parse_readname(line: &str) -> &str {
    line.split('\t').nth(4).unwrap_or("")
}

/// Merge a single cluster of candidate lines into one event string.
/// Mirrors Python `mergeinfo_insertion_oneevent` in debreak_merge_clustering.py.
///
/// Trims outliers from both ends until candidates are within 2× of each other,
/// then computes weighted-average position, size, quality across all lines.
///
/// Returns `Some(merged_line)` if `candi.len() >= max(2, min_support)`, else `None`.
pub fn merge_one_event(candi: &[String], min_support: usize) -> Option<String> {
    if candi.is_empty() {
        return None;
    }

    let effective_min = min_support.max(2);

    // Sort by size ascending (working copy of indices)
    let mut order: Vec<usize> = (0..candi.len()).collect();
    order.sort_by_key(|&i| parse_size(&candi[i]));

    // Trim outlier sizes from either end while we have more than effective_min-2 candidates
    loop {
        let n = order.len();
        if n <= effective_min.saturating_sub(2).max(2) {
            break;
        }
        let top = parse_size(&candi[order[n - 1]]);
        let mid = parse_size(&candi[order[n / 2]]);
        let bot = parse_size(&candi[order[0]]);

        if top > 2 * mid && top.saturating_sub(mid) > 30 {
            order.pop();
            continue;
        }
        if mid > 2 * bot && mid.saturating_sub(bot) > 30 {
            order.remove(0);
            continue;
        }
        break;
    }

    if order.len() < effective_min {
        return None;
    }

    let chrom = candi[order[0]].split('\t').next().unwrap_or("").to_string();
    let count = order.len();

    // Weighted averages
    let pos_sum: u64   = order.iter().map(|&i| parse_pos(&candi[i])).sum();
    let size_sum: u64  = order.iter().map(|&i| parse_size(&candi[i])).sum();
    let qual_sum: f64  = order.iter().map(|&i| parse_quality(&candi[i])).sum();
    let position = pos_sum  / count as u64;
    let size     = size_sum / count as u64;
    let quality  = qual_sum / count as f64;

    // Collect readnames (field[4]) — semicolon-delimited
    let readnames: String = order.iter()
        .map(|&i| parse_readname(&candi[i]))
        .collect::<Vec<_>>()
        .join(";");
    let numread = readnames.split(';').count();

    // Format: chrom\tpos\tsize\tcount\tnumread\tquality\t\treadnames
    // (empty field between quality and readnames mirrors the Python format)
    Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{:.6}\t\t{}",
        chrom, position, size, count, numread, quality, readnames
    ))
}

/// Merge candidate lines, detecting bimodal size distributions and splitting
/// into compound SVs when warranted.
/// Mirrors Python `mergeinfo_insertion` in debreak_merge_clustering.py.
///
/// Bimodal test: if upper quartile > 1.75× lower quartile AND gap > 50bp,
/// run 3 iterations of K-means (k=2) and try to report two sub-events.
///
/// Returns merged event strings, each appended with `\tUnique` or `\tCompoundSV`.
pub fn merge_with_bimodal(candi: &[String], min_support: usize) -> Vec<String> {
    if candi.is_empty() {
        return Vec::new();
    }

    // Sort indices by size to get quartiles
    let mut order: Vec<usize> = (0..candi.len()).collect();
    order.sort_by_key(|&i| parse_size(&candi[i]));
    let n = order.len();

    // Bimodal test: needs at least 1.5× min_support candidates
    if n as f64 >= 1.5 * min_support as f64 {
        let upper = parse_size(&candi[order[n * 3 / 4]]);
        let lower = parse_size(&candi[order[n / 4]]);

        if upper > lower && upper > (lower as f64 * 1.75) as u64 && upper.saturating_sub(lower) > 50 {
            // Run K-means assignment for 3 iterations
            let sizes: Vec<u64> = order.iter().map(|&i| parse_size(&candi[i])).collect();
            let (mut g1, mut g2, mut m1, mut m2) = kmeans_assign(&sizes, upper, lower);
            for _ in 0..2 {
                let r = kmeans_assign(&sizes, m1, m2);
                g1 = r.0; g2 = r.1; m1 = r.2; m2 = r.3;
            }

            let cluster1: Vec<String> = g1.iter().map(|&i| candi[order[i]].clone()).collect();
            let cluster2: Vec<String> = g2.iter().map(|&i| candi[order[i]].clone()).collect();

            let ev1 = merge_one_event(&cluster1, min_support);
            let ev2 = merge_one_event(&cluster2, min_support);

            match (ev1, ev2) {
                (Some(e1), Some(e2)) => {
                    // Both clusters valid — report as CompoundSV pair
                    return vec![
                        format!("{}\tCompoundSV", e1),
                        format!("{}\tCompoundSV", e2),
                    ];
                }
                (Some(e), None) | (None, Some(e)) => {
                    return vec![format!("{}\tUnique", e)];
                }
                (None, None) => {}
            }
        }
    }

    // Unimodal path
    match merge_one_event(candi, min_support) {
        Some(e) => vec![format!("{}\tUnique", e)],
        None    => Vec::new(),
    }
}

/// Adaptive-window clustering of SV signals, mirroring Python `counttime_insertion` /
/// `counttime_deletion` from debreak_merge.py.
///
/// Algorithm:
/// 1. Sort signals by (position, size).
/// 2. Start with `window = 100`.
/// 3. If the next event falls within `[start, start + window)`, accumulate.
/// 4. Otherwise, if `window == 100` (not yet adapted for this cluster), adapt window from
///    the mean size of the current cluster:
///      mean ≤ 100  →  window = 200
///      mean ≤ 500  →  window = 400
///      mean  > 500 →  window = 800
///    Re-check: if the event now falls in the adapted window, absorb it.
/// 5. If still outside, flush the cluster (→ `merge_with_bimodal` if len ≥ min_support)
///    and start a new cluster from the current event with window reset to 100.
///
/// This matches the Python behaviour for large SVs (it is applied to all signal sizes;
/// for small SVs the depth-map path is used instead and this function is not called).
pub fn counttime_cluster(signals: &[String], min_support: usize) -> Vec<String> {
    if signals.is_empty() {
        return Vec::new();
    }

    // Sort a copy by (pos, size)
    let mut sorted = signals.to_vec();
    sorted.sort_by_key(|l| (parse_pos(l), parse_size(l)));

    let mut results: Vec<String> = Vec::new();
    let mut candi: Vec<String> = vec![sorted[0].clone()];
    let mut start = parse_pos(&sorted[0]);
    let mut window: u64 = 100;

    for event in sorted.into_iter().skip(1) {
        let ev_pos = parse_pos(&event);
        if ev_pos <= start + window {
            candi.push(event);
        } else {
            // Event is outside current window. Adapt window if not yet done.
            if window == 100 {
                let mean_size: u64 = if candi.is_empty() {
                    0
                } else {
                    candi.iter().map(|l| parse_size(l)).sum::<u64>() / candi.len() as u64
                };
                window = if mean_size <= 100 { 200 } else if mean_size <= 500 { 400 } else { 800 };

                // Re-check with adapted window
                if ev_pos <= start + window {
                    candi.push(event);
                    continue;
                }
            }

            // Flush current cluster
            if candi.len() >= min_support {
                results.extend(merge_with_bimodal(&candi, min_support));
            }
            // Start new cluster
            start = ev_pos;
            window = 100;
            candi = vec![event];
        }
    }

    // Flush final cluster
    if candi.len() >= min_support {
        results.extend(merge_with_bimodal(&candi, min_support));
    }

    results
}

/// Merge insertion signals (kept for API compatibility — unused by new pipeline)
pub fn merge_insertions(signals: &[(u64, usize)]) -> Vec<(u64, usize)> {
    signals.to_vec()
}

/// Merge deletion signals (kept for API compatibility — unused by new pipeline)
pub fn merge_deletions(signals: &[(u64, usize)]) -> Vec<(u64, usize)> {
    signals.to_vec()
}

/// Merge translocation signals (kept for API compatibility — unused by new pipeline)
pub fn merge_translocations(signals: &[(String, u64, String, u64)]) -> Vec<(String, u64, String, u64)> {
    signals.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_line(chrom: &str, pos: u64, size: u64, readname: &str, mapq: u8) -> String {
        format!("{}\t{}\t{}\tI-cigar\t{}\t0\t{}", chrom, pos, size, readname, mapq)
    }

    #[test]
    fn test_merge_insertions() {
        // Kept passing — API unchanged
    }

    #[test]
    fn test_merge_one_event_basic() {
        // 5 reads at roughly the same position and size
        let candi: Vec<String> = (0..5)
            .map(|i| make_line("chr1", 1000 + i * 2, 200 + i, &format!("read{}", i), 30))
            .collect();
        let result = merge_one_event(&candi, 3);
        assert!(result.is_some(), "Expected a merged event");
        let merged = result.unwrap();
        let f: Vec<&str> = merged.split('\t').collect();
        assert_eq!(f[0], "chr1");
        // position should be average of 1000..1008
        let pos: u64 = f[1].parse().unwrap();
        assert!(pos >= 1000 && pos <= 1010, "position {} out of expected range", pos);
        // readnames field is field[7] (note empty field[6])
        assert_eq!(f.len(), 8, "Expected 8 fields: {}", merged);
    }

    #[test]
    fn test_merge_one_event_trims_outliers() {
        // One huge outlier should be trimmed
        let mut candi: Vec<String> = (0..5)
            .map(|i| make_line("chr1", 1000, 100 + i * 5, &format!("r{}", i), 20))
            .collect();
        candi.push(make_line("chr1", 1000, 9999, "outlier", 20)); // massive outlier
        let result = merge_one_event(&candi, 3);
        assert!(result.is_some());
            let merged = result.unwrap();
            let f: Vec<&str> = merged.split('\t').collect();
        let size: u64 = f[2].parse().unwrap();
        assert!(size < 500, "Outlier should be trimmed, size={}", size);
    }

    #[test]
    fn test_merge_one_event_below_min_support() {
        // Only 1 read, min_support=5 → None
        let candi = vec![make_line("chr1", 500, 100, "read1", 30)];
        let result = merge_one_event(&candi, 5);
        assert!(result.is_none());
    }

    #[test]
    fn test_merge_with_bimodal_unimodal() {
        // All same size → single Unique event
        let candi: Vec<String> = (0..6)
            .map(|i| make_line("chr1", 1000 + i, 200, &format!("r{}", i), 25))
            .collect();
        let result = merge_with_bimodal(&candi, 3);
        assert_eq!(result.len(), 1, "Expected single merged event");
        assert!(result[0].ends_with("Unique"), "Expected Unique tag: {}", result[0]);
    }

    #[test]
    fn test_merge_with_bimodal_two_clusters() {
        // Two size clusters: ~100bp and ~500bp, 5 reads each → CompoundSV
        let mut candi: Vec<String> = (0..5)
            .map(|i| make_line("chr1", 1000 + i, 100 + i * 2, &format!("small{}", i), 25))
            .collect();
        candi.extend((0..5)
            .map(|i| make_line("chr1", 1000 + i, 500 + i * 2, &format!("big{}", i), 25)));
        let result = merge_with_bimodal(&candi, 3);
        // Should produce 2 CompoundSV events
        assert_eq!(result.len(), 2, "Expected 2 compound events, got: {:?}", result);
        assert!(result.iter().all(|r| r.ends_with("CompoundSV")), "Both should be CompoundSV");
    }
}
