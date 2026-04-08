/// Structural variant clustering, genotyping, and reconciliation
/// Ported from Python debreak_merge_clustering.py (cluster, cluster_ins, genotype, filterae)

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use anyhow::Result;
use log::info;
use crate::merge::merge_ops::merge_with_bimodal;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn read_debreak_temp(outpath: &str, chrom: &str) -> Result<Vec<String>> {
    let path = format!("{}debreak_workspace/read_to_contig_{}.debreak.temp", outpath, chrom);
    match File::open(&path) {
        Err(_) => Ok(Vec::new()),
        Ok(f) => Ok(BufReader::new(f)
            .lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.is_empty())
            .collect()),
    }
}

fn sig_pos(line: &str) -> u64 {
    line.split('\t').nth(1).and_then(|s| s.parse().ok()).unwrap_or(0)
}

fn sig_size(line: &str) -> u64 {
    line.split('\t').nth(2).and_then(|s| s.parse().ok()).unwrap_or(0)
}

/// Fixed-window grouping for large SVs (>2000 bp).
fn cluster_large_svs(signals: &[String], min_support: usize, window: u64) -> Vec<String> {
    if signals.is_empty() { return Vec::new(); }
    let mut sorted = signals.to_vec();
    sorted.sort_by_key(|l| (sig_pos(l), sig_size(l)));
    let mut results = Vec::new();
    let mut candi: Vec<String> = Vec::new();
    let mut start = sig_pos(&sorted[0]);
    for event in sorted {
        if sig_pos(&event) <= start + window {
            candi.push(event);
        } else {
            if candi.len() >= min_support {
                results.extend(merge_with_bimodal(&candi, min_support));
            }
            start = sig_pos(&event);
            candi = vec![event];
        }
    }
    if candi.len() >= min_support {
        results.extend(merge_with_bimodal(&candi, min_support));
    }
    results
}

/// Depth-map spatial clustering for small SVs (<=3000bp).
/// mark_fn(pos,size) -> (start,end) marks the depth array.
/// assign_fn(pos,size,r_start,r_end) -> bool assigns signal to region.
/// Returns merged lines with "\t\tr_start\tr_end" appended (stripped by caller).
fn cluster_small_svs_depth_map(
    signals: &[String],
    contig_length: u64,
    max_depth: usize,
    min_support: usize,
    mark_fn: impl Fn(u64, u64) -> (usize, usize),
    assign_fn: impl Fn(u64, u64, usize, usize) -> bool,
) -> Vec<String> {
    if signals.is_empty() || contig_length == 0 { return Vec::new(); }
    let len = contig_length as usize;
    let mut depth: Vec<u32> = vec![0u32; len];
    for sig in signals {
        let (s, e) = mark_fn(sig_pos(sig), sig_size(sig));
        let s = s.min(len.saturating_sub(1));
        let e = e.min(len);
        if s < e { for v in &mut depth[s..e] { *v += 1; } }
    }
    let threshold = 3u32;
    let mut regions: Vec<(usize, usize)> = Vec::new();
    let mut in_block = false;
    let mut block_start = 0usize;
    let mut max_dep = 0u32;
    for i in 0..len {
        let d = depth[i];
        if in_block {
            let floor = (max_dep as f64 / 10.0).max(threshold as f64) as u32;
            if d >= floor {
                if d > max_dep { max_dep = d; }
            } else {
                if (max_dep as usize) <= max_depth { regions.push((block_start, i)); }
                in_block = false; max_dep = 0;
            }
        } else if d > threshold {
            in_block = true; block_start = i; max_dep = d;
        }
    }
    if in_block && (max_dep as usize) <= max_depth { regions.push((block_start, len)); }
    let mut sv_results: Vec<String> = Vec::new();
    for (r_start, r_end) in &regions {
        let region_sigs: Vec<String> = signals.iter()
            .filter(|sig| assign_fn(sig_pos(sig), sig_size(sig), *r_start, *r_end))
            .cloned().collect();
        for m in merge_with_bimodal(&region_sigs, min_support) {
            sv_results.push(format!("{}\t\t{}\t{}", m, r_start, r_end));
        }
    }
    sv_results
}

fn strip_region_suffix(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() > 3 { fields[..fields.len() - 3].join("\t") } else { line.to_string() }
}

fn dedup_large_small(large_merged: Vec<String>, small_merged: &[String], ins_mode: bool) -> Vec<String> {
    let mut kept = Vec::new();
    for large in large_merged {
        let l_pos = sig_pos(&large); let l_size = sig_size(&large); let l_end = l_pos + l_size;
        let overlaps = small_merged.iter().any(|s| {
            let s_pos = sig_pos(s); let s_size = sig_size(s);
            let s_end = if ins_mode { s_pos + 1 } else { s_pos + s_size };
            l_end.min(s_end) > l_pos.max(s_pos)
                && (s_size as f64 * 0.8) <= l_size as f64
                && l_size as f64 <= (s_size as f64 / 0.8)
        });
        if !overlaps { kept.push(large); }
    }
    kept
}

// ---------------------------------------------------------------------------
// Public clustering functions
// ---------------------------------------------------------------------------

/// Cluster deletion signals → del_merged_{chrom}.
/// Mirrors Python debreak_merge_clustering.cluster().
pub fn cluster(
    outpath: &str,
    chrom: &str,
    contig_length: u64,
    min_support: usize,
    max_depth: usize,
) -> Result<()> {
    info!("Clustering deletions for {}", chrom);
    let all_sv = read_debreak_temp(outpath, chrom)?;

    let large_del: Vec<String> = all_sv.iter().filter(|l| l.contains("D-") && sig_size(l) > 2000).cloned().collect();
    let large_merged = cluster_large_svs(&large_del, min_support, 1600);

    let small_del: Vec<String> = all_sv.iter().filter(|l| l.contains("D-") && sig_size(l) <= 3000).cloned().collect();
    let small_annotated = cluster_small_svs_depth_map(&small_del, contig_length, max_depth, min_support,
        |pos, size| { let s = (pos as usize).saturating_sub(1); let e = (pos + size) as usize; (s, e) },
        |pos, size, r_start, r_end| { let end = pos + size; end > r_start as u64 && pos < r_end as u64 },
    );
    let small_merged: Vec<String> = small_annotated.iter().map(|l| strip_region_suffix(l)).collect();

    let mut final_svs = dedup_large_small(large_merged, &small_merged, false);
    final_svs.extend(small_merged);
    final_svs.sort_by_key(|l| (sig_pos(l), sig_size(l)));

    let output_file = format!("{}ae_merge_workspace/del_merged_{}", outpath, chrom);
    let mut output = File::create(&output_file)?;
    for sv in &final_svs { writeln!(output, "{}", sv)?; }
    info!("del cluster {}: {} signals -> {} events", chrom, all_sv.len(), final_svs.len());
    Ok(())
}

/// Cluster insertion/inversion signals → ins_merged_{chrom} or inv_merged_{chrom}.
/// Mirrors Python debreak_merge_clustering.cluster_ins().
pub fn cluster_insertions(
    outpath: &str,
    chrom: &str,
    contig_length: u64,
    min_support: usize,
    max_depth: usize,
    sv_type: &str,
) -> Result<()> {
    let signal_marker = if sv_type == "ins" { "I-" } else { "INV-" };
    let type_label    = if sv_type == "ins" { "insertion" } else { "inversion" };
    info!("Clustering {}s for {}", type_label, chrom);
    let all_sv = read_debreak_temp(outpath, chrom)?;

    let large_sv: Vec<String> = all_sv.iter().filter(|l| l.contains(signal_marker) && sig_size(l) > 2000).cloned().collect();
    let large_merged = cluster_large_svs(&large_sv, min_support, 1600);

    let small_sv: Vec<String> = all_sv.iter().filter(|l| l.contains(signal_marker) && sig_size(l) <= 3000).cloned().collect();
    let small_annotated = cluster_small_svs_depth_map(&small_sv, contig_length, max_depth, min_support,
        |pos, _size| { let s = (pos as usize).saturating_sub(101); let e = pos as usize + 101; (s, e) },
        |pos, _size, r_start, r_end| {
            let sig_s = (pos as usize).saturating_sub(50);
            let sig_e = pos as usize + 50;
            sig_e > r_start && sig_s < r_end
        },
    );
    let small_merged: Vec<String> = small_annotated.iter().map(|l| strip_region_suffix(l)).collect();

    let mut final_svs = dedup_large_small(large_merged, &small_merged, true);
    final_svs.extend(small_merged);
    final_svs.sort_by_key(|l| (sig_pos(l), sig_size(l)));

    let prefix = if sv_type == "ins" { "ins" } else { "inv" };
    let output_file = format!("{}ae_merge_workspace/{}_merged_{}", outpath, prefix, chrom);
    let mut output = File::create(&output_file)?;
    for sv in &final_svs { writeln!(output, "{}", sv)?; }
    info!("{} cluster {}: {} signals -> {} events", type_label, chrom, all_sv.len(), final_svs.len());
    Ok(())
}

// ---------------------------------------------------------------------------
// Genotyping
// ---------------------------------------------------------------------------

/// Count reads overlapping [start, end) via samtools view -c.
fn count_reads_in_region(bam: &str, chrom: &str, start: u64, end: u64) -> u64 {
    if end <= start {
        return 0;
    }
    let region = format!("{}:{}-{}", chrom, start + 1, end);
    match std::process::Command::new("samtools")
        .args(["view", "-c", bam, &region])
        .output()
    {
        Ok(o) => String::from_utf8_lossy(&o.stdout)
            .trim()
            .parse::<u64>()
            .unwrap_or(0),
        Err(_) => 0,
    }
}

/// Assign 1/1 or 1/0 genotypes via flanking depth. Writes assembly_errors.bed-gt.
/// Mirrors Python debreak_merge_clustering.genotype().
pub fn genotype(coverage: usize, outpath: &str) -> Result<()> {
    info!("Genotyping structural errors (avg coverage: {})", coverage);
    let bed_path = format!("{}assembly_errors.bed", outpath);
    let gt_path  = format!("{}assembly_errors.bed-gt", outpath);
    let bed_file = match File::open(&bed_path) {
        Err(_) => { File::create(&gt_path)?; return Ok(()); }
        Ok(f) => f,
    };
    let bam = format!("{}read_to_contig.bam", outpath);
    let mut gt_out = File::create(&gt_path)?;
    for line in BufReader::new(bed_file).lines() {
        let line = line?;
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 { continue; }
        let chrom = fields[0];
        let start: i64 = fields[1].parse().unwrap_or(0);
        let stop:  i64 = fields[2].parse().unwrap_or(0);
        let supporting: f64 = fields[3].parse().unwrap_or(0.0);
        if start < 0 { continue; }
        let left_cov  = count_reads_in_region(&bam, chrom, (start - 100).max(0) as u64, start as u64);
        let right_cov = count_reads_in_region(&bam, chrom, stop as u64, (stop + 100) as u64);
        let min_cov = left_cov.min(right_cov);
        let gt = if supporting >= 0.6 * min_cov as f64 { "1/1" } else { "1/0" };
        // Match Python output columns: original fields + gt + left + right + min
        writeln!(gt_out, "{}\t{}\t{}\t{}\t{}", line, gt, left_cov, right_cov, min_cov)?;
    }
    info!("Genotyping complete");
    Ok(())
}

// ---------------------------------------------------------------------------
// Reconciliation stub (Todo 7)
// ---------------------------------------------------------------------------

/// Filter and reconcile expansion/collapse pairs.
/// Mirrors Python filterae() in debreak_merge_clustering.py.
pub fn filter_errors(coverage: usize, outpath: &str, min_size: usize, datatype: &str) -> Result<usize> {
    info!("Filtering errors (avg coverage: {})", coverage);

    let gt_path = format!("{}assembly_errors.bed-gt", outpath);
    let allsv_raw: Vec<String> = match File::open(&gt_path) {
        Ok(f) => BufReader::new(f)
            .lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.is_empty())
            .collect(),
        Err(_) => {
            File::create(format!("{}structural_error.bed", outpath))?;
            return Ok(0);
        }
    };

    let rat = if datatype == "hifi" { 0.8 } else { 0.7 };
    let highcov = (coverage * 2) as i64;
    let lowcov = (coverage / 2) as i64;

    let mut exp: Vec<String> = allsv_raw.iter().filter(|c| c.contains("Exp")).cloned().collect();
    let mut col: Vec<String> = allsv_raw.iter().filter(|c| c.contains("Col")).cloned().collect();
    let inv: Vec<String> = allsv_raw.iter().filter(|c| c.contains("Inv")).cloned().collect();

    let uniq_count = |reads: &str| -> usize {
        let mut v: Vec<&str> = reads.split(';').filter(|s| !s.is_empty()).collect();
        v.sort_unstable();
        v.dedup();
        v.len()
    };
    let uniq_join = |reads: &str| -> String {
        let mut v: Vec<&str> = reads.split(';').filter(|s| !s.is_empty()).collect();
        v.sort_unstable();
        v.dedup();
        v.join(";")
    };
    let parse_size_field = |s: &str| -> Vec<i64> {
        let raw = s.strip_prefix("Size=").unwrap_or(s);
        raw.split(';').filter_map(|x| x.parse::<i64>().ok()).collect()
    };

    let mut reconciled: Vec<String> = Vec::new();
    let mut exp_only: Vec<String> = Vec::new();

    for e in exp.drain(..) {
        let c: Vec<&str> = e.split('\t').collect();
        if c.len() < 11 { continue; }
        let chrom = c[0];
        let c_start: i64 = c[1].parse().unwrap_or(0);
        let c_end: i64 = c[2].parse().unwrap_or(0);
        let c_sup: i64 = c[3].parse().unwrap_or(0);
        let c_size = parse_size_field(c[5]).get(0).copied().unwrap_or(0);

        let mut matched = false;
        let mut remove_idx: Option<usize> = None;

        for (j, dline) in col.iter().enumerate() {
            let d: Vec<&str> = dline.split('\t').collect();
            if d.len() < 11 { continue; }
            if d[0] != chrom { continue; }

            let d_start: i64 = d[1].parse().unwrap_or(0);
            let in_window = c_start - 250 <= d_start && d_start <= c_end + 250;
            let d_size = parse_size_field(d[5]).get(0).copied().unwrap_or(0);
            if !in_window || c_size >= 20 * d_size { continue; }

            matched = true;
            let expreads = c[6];
            let colreads = d[6];
            let goodexp = uniq_count(expreads) as i64;
            let goodcol = uniq_count(colreads) as i64;
            let totaln = {
                let mut both: Vec<&str> = expreads.split(';').chain(colreads.split(';')).filter(|s| !s.is_empty()).collect();
                both.sort_unstable();
                both.dedup();
                both.len() as i64
            };

            let d_sup: i64 = d[3].parse().unwrap_or(1).max(1);
            let ratio = c_sup as f64 / d_sup as f64;

            if (0.33..=3.0).contains(&ratio) {
                let thresh = (goodexp + goodcol / 2).min(goodcol + goodexp / 2);
                if totaln < thresh {
                    if c_size > d_size + min_size as i64 {
                        reconciled.push(format!(
                            "{}\t{}\t{}\t{}\tExpansion\tSize={}\t{}\t{}\t{}\t{}\t{}",
                            chrom, c[1], c[2], goodexp, c_size - d_size, c[7], c[8], c[9], c[10], uniq_join(expreads)
                        ));
                    }
                    if c_size < d_size - min_size as i64 {
                        reconciled.push(format!(
                            "{}\t{}\t{}\t{}\tCollapse\tSize={}\t{}\t{}\t{}\t{}\t{}",
                            d[0], d[1], d[2], goodcol, d_size - c_size, d[7], d[8], d[9], d[10], uniq_join(colreads)
                        ));
                    }
                } else {
                    reconciled.push(format!(
                        "{}\t{};{}\t{};{}\t{}\tHaplotypeSwitch\tSize={};{}\t-/-\t{}\t{}\t{}\t{}:{}\t{};{}",
                        chrom, c[1], d[1], c[2], d[2], totaln, c_end - c_start, d_size,
                        c[8], c[9], c[10], uniq_join(expreads), uniq_join(colreads), goodexp, goodcol
                    ));
                }
            } else if ratio < 0.33 {
                reconciled.push(format!(
                    "{}\t{}\t{}\t{}\tCollapse\tSize={};{}\t-/-\t{}\t{}\t{}\t{}",
                    d[0], d[1], d[2], goodcol, c_end - c_start, d_size, d[8], d[9], d[10], uniq_join(colreads)
                ));
            } else {
                reconciled.push(format!(
                    "{}\t{}\t{}\t{}\tExpansion\tSize={};{}\t-/-\t{}\t{}\t{}\t{}",
                    chrom, c[1], c[2], goodexp, c_end - c_start, d_size, c[8], c[9], c[10], uniq_join(expreads)
                ));
            }

            remove_idx = Some(j);
            break;
        }

        if let Some(j) = remove_idx { col.remove(j); }
        if !matched { exp_only.push(e); }
    }

    let mut allsv: Vec<String> = reconciled;

    for line in exp_only.into_iter().chain(col.into_iter()) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 { continue; }
        let reads = uniq_join(f[6]);
        let good = uniq_count(f[6]);
        allsv.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f[0], f[1], f[2], good, f[4], f[5], f[7], f[8], f[9], f[10], reads
        ));
    }

    for line in inv {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 { continue; }
        let reads = uniq_join(f[5]);
        let good = uniq_count(f[5]);
        let inv_size: i64 = f[2].parse::<i64>().unwrap_or(0) - f[1].parse::<i64>().unwrap_or(0);
        allsv.push(format!(
            "{}\t{}\t{}\t{}\t{}\tSize={}\t{}\t{}\t{}\t{}\t{}",
            f[0], f[1], f[2], good, f[4], inv_size.max(0), f[6], f[7], f[8], f[9], reads
        ));
    }

    let mut filtered: Vec<String> = Vec::new();
    for line in allsv {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 { continue; }
        let support: i64 = f[3].parse().unwrap_or(0);
        let depth_min: i64 = f[9].parse().unwrap_or(0);
        let sizes = parse_size_field(f[5]);
        let max_size = sizes.into_iter().max().unwrap_or(0);
        if max_size < min_size as i64 { continue; }
        if support >= 10 && (support as f64) >= rat * depth_min as f64 && lowcov <= depth_min && depth_min < highcov {
            filtered.push(line);
        }
    }

    let structural_path = format!("{}structural_error.bed", outpath);
    let mut out = File::create(&structural_path)?;
    writeln!(out, "#Contig_Name\tStart_Position\tEnd_Position\tSupporting_Read\tType\tSize\tHaplotype_Info\tDepth_Left\tDepth_Right\tDepth_Min\tSupporting_Read_Name\tHaplotype_Switch_Info")?;
    for line in &filtered { writeln!(out, "{}", line)?; }

    let exp_n = filtered.iter().filter(|l| l.contains("Exp")).count();
    let col_n = filtered.iter().filter(|l| l.contains("Col")).count();
    let het_n = filtered.iter().filter(|l| l.contains("Haplo")).count();
    let inv_n = filtered.iter().filter(|l| l.contains("Inv")).count();
    let mut summary = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(format!("{}summary_statistics.txt", outpath))?;
    writeln!(summary, "Structural error\t{}", filtered.len())?;
    writeln!(summary, "Expansion\t{}", exp_n)?;
    writeln!(summary, "Collapse\t{}", col_n)?;
    writeln!(summary, "Haplotype switch\t{}", het_n)?;
    writeln!(summary, "Inversion\t{}", inv_n)?;

    let _ = std::fs::remove_file(format!("{}assembly_errors.bed", outpath));
    let _ = std::fs::remove_file(format!("{}assembly_errors.bed-gt", outpath));
    let _ = std::fs::remove_file(format!("{}read_to_contig.debreak.temp", outpath));

    let mut totalbase = 0usize;
    for line in &filtered {
        if line.contains("Inv") {
            totalbase += 200;
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 { continue; }
        let sizes = parse_size_field(f[5]);
        let size = sizes.into_iter().min().unwrap_or(0).min(10000);
        totalbase += size.max(0) as usize;
    }

    info!("Total filtered errors: {}", filtered.len());
    Ok(totalbase)
}

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------

/// Write structural-error-inspector.tsv from structural_error.bed (or fallback to ae_merge_workspace).
pub fn write_structural_error_tsv(outpath: &str, base_name: &str) -> Result<()> {
    info!("Writing structural error inspector TSV");
    let bed_path = format!("{}structural_error.bed", outpath);
    let lines: Vec<String> = if let Ok(f) = File::open(&bed_path) {
        BufReader::new(f).lines().filter_map(|l| l.ok())
            .filter(|l| !l.is_empty() && !l.starts_with('#')).collect()
    } else {
        let mut all = Vec::new();
        if let Ok(entries) = std::fs::read_dir(format!("{}ae_merge_workspace/", outpath)) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(n) = path.file_name().and_then(|n| n.to_str()) {
                    if n.starts_with("del_merged_") || n.starts_with("ins_merged_") || n.starts_with("inv_merged_") {
                        if let Ok(f) = File::open(&path) {
                            all.extend(BufReader::new(f).lines().filter_map(|l| l.ok()).filter(|l| !l.is_empty()));
                        }
                    }
                }
            }
        }
        all.sort_by(|a, b| {
            let pa = (a.split('\t').next().unwrap_or(""), sig_pos(a));
            let pb = (b.split('\t').next().unwrap_or(""), sig_pos(b));
            pa.cmp(&pb)
        });
        all
    };
    let tsv_path = format!("{}{}-structural-error-inspector.tsv", outpath, base_name);
    let mut tsv = File::create(&tsv_path)?;
    writeln!(tsv, "#Contig_Name\tStart_Position\tEnd_Position\tSupporting_Read\tType\tSize\tHaplotype_Info\tDepth_Left\tDepth_Right\tDepth_Min\tSupporting_Read_Name\tHaplotype_Switch_Info")?;
    for sv in &lines { writeln!(tsv, "{}", sv)?; }
    info!("Wrote {} records to {}", lines.len(), tsv_path);
    Ok(())
}

/// Append extended summary statistics to summary_statistics file.
pub fn write_summary_statistics_extended(
    outpath: &str,
    fasta_stats: &crate::static_analysis::FastaStats,
    coverage: usize,
    ae_len_structural_error: usize,
    ae_len_base_error: usize,
) -> Result<()> {
    let total_errors = ae_len_structural_error + ae_len_base_error;
    let valid_bases  = fasta_stats.total_length;
    let long_read_qv = if total_errors == 0 || valid_bases == 0 { f64::INFINITY }
        else { -10.0 * (total_errors as f64 / valid_bases as f64).log10() };
    let mut summary = std::fs::OpenOptions::new().create(true).append(true)
        .open(format!("{}summary_statistics", outpath))?;
    writeln!(summary)?;
    writeln!(summary, "Assembly error size\t{}", total_errors)?;
    writeln!(summary, "Structural error size\t{}", ae_len_structural_error)?;
    writeln!(summary, "Base error\t{}", ae_len_base_error)?;
    writeln!(summary, "Average coverage\t{}", coverage)?;
    if long_read_qv.is_infinite() { writeln!(summary, "long_read_QV\tInf")?; }
    else { writeln!(summary, "long_read_QV\t{:.2}", long_read_qv)?; }
    writeln!(summary)?;
    writeln!(summary, "Total contigs\t{}", fasta_stats.chromosomes.len())?;
    writeln!(summary, "Large contigs (>=10Mbp)\t{}", fasta_stats.chromosomes_large.len())?;
    writeln!(summary, "Total length\t{}", fasta_stats.total_length)?;
    writeln!(summary, "Longest contig\t{}", fasta_stats.largest_contig_length)?;
    writeln!(summary, "N50\t{}", fasta_stats.n50)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_large_svs_groups_correctly() {
        let lines: Vec<String> = vec![
            "chr1\t1000\t3000\tD-cigar\tread1\t0\t30".to_string(),
            "chr1\t1500\t3000\tD-cigar\tread2\t0\t30".to_string(),
            "chr1\t2800\t3000\tD-cigar\tread3\t0\t30".to_string(),
            "chr1\t1200\t3000\tD-cigar\tread4\t0\t30".to_string(),
            "chr1\t6000\t3000\tD-cigar\tread5\t0\t30".to_string(),
        ];
        let result = cluster_large_svs(&lines, 3, 1600);
        assert_eq!(result.len(), 1, "Expected 1 merged group, got {:?}", result);
    }

    #[test]
    fn test_cluster_large_svs_below_support() {
        let lines: Vec<String> = vec![
            "chr1\t1000\t3000\tD-cigar\tread1\t0\t30".to_string(),
            "chr1\t1100\t3000\tD-cigar\tread2\t0\t30".to_string(),
        ];
        let result = cluster_large_svs(&lines, 5, 1600);
        assert!(result.is_empty());
    }

    #[test]
    fn test_strip_region_suffix() {
        let line = "chr1\t1000\t200\t5\t5\t30.0\tUnique\t\t100\t200";
        let stripped = strip_region_suffix(line);
        assert!(!stripped.ends_with("100\t200"), "Region suffix should be removed: {}", stripped);
        assert!(stripped.contains("chr1"), "Chrom should remain");
    }

    #[test]
    fn test_cluster() {}
}
