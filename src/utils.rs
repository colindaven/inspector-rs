/// Utility functions used across modules

pub mod bam;

use std::path::Path;
use std::fs;
use anyhow::Result;

/// Check if file exists and is readable
pub fn validate_input_file(path: &str) -> Result<()> {
    let p = Path::new(path);
    if !p.exists() {
        anyhow::bail!("Input file not found: {}", path);
    }
    if p.is_dir() {
        anyhow::bail!("Expected a file but found a directory: {}", path);
    }
    // Attempt to open for reading to catch permission errors on read-only filesystems
    fs::File::open(p).map_err(|e| anyhow::anyhow!("Cannot read file '{}': {}", path, e))?;
    Ok(())
}

/// Check if all input files in a list exist
pub fn validate_input_files(paths: &[String]) -> Result<()> {
    for path in paths {
        validate_input_file(path)?;
    }
    Ok(())
}

/// Ensure output directory exists and is writable
pub fn ensure_output_dir(path: &str) -> Result<()> {
    if !Path::new(path).exists() {
        fs::create_dir_all(path)
            .map_err(|e| anyhow::anyhow!("Cannot create output directory '{}': {}", path, e))?;
    }
    // Verify we can write to it
    let test_file = format!("{}.inspector_write_test", path.trim_end_matches('/'));
    match fs::write(&test_file, b"") {
        Ok(_) => { let _ = fs::remove_file(&test_file); }
        Err(e) => anyhow::bail!("Output directory is not writable '{}': {}", path, e),
    }
    Ok(())
}

/// Check if an external tool is available in PATH
pub fn require_tool(tool: &str) -> Result<String> {
    let output = std::process::Command::new(tool)
        .arg("--version")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output()
        .map_err(|_| anyhow::anyhow!(
            "Required tool '{}' not found in PATH. Please install it and ensure it is on your PATH.",
            tool
        ))?;
    // Version strings are often on stdout; some tools (e.g. samtools) print to stderr
    let version_line = String::from_utf8_lossy(&output.stdout)
        .lines()
        .next()
        .unwrap_or("")
        .trim()
        .to_string();
    let version_line = if version_line.is_empty() {
        String::from_utf8_lossy(&output.stderr)
            .lines()
            .next()
            .unwrap_or("unknown version")
            .trim()
            .to_string()
    } else {
        version_line
    };
    Ok(version_line)
}

/// Custom writer that writes to both stdout and a file
struct TeeWriter {
    file: std::sync::Arc<std::sync::Mutex<fs::File>>,
}

impl TeeWriter {
    fn new(file: fs::File) -> Self {
        TeeWriter {
            file: std::sync::Arc::new(std::sync::Mutex::new(file)),
        }
    }
}

impl std::io::Write for TeeWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        use std::io::Write;
        
        // Write to stdout
        std::io::stdout().write_all(buf)?;
        std::io::stdout().flush()?;
        
        // Write to file
        if let Ok(mut f) = self.file.lock() {
            f.write_all(buf)?;
            f.flush()?;
        }
        
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        std::io::stdout().flush()?;
        if let Ok(mut f) = self.file.lock() {
            f.flush()?;
        }
        Ok(())
    }
}

/// Set up file logger to write all messages to Inspector.log AND stdout
/// If truncate=true, overwrites existing log file; if false, appends
pub fn init_file_logger(log_file_path: &str, truncate: bool) -> Result<()> {
    use log::LevelFilter;
    use std::io::Write;

    let log_path = Path::new(log_file_path);
    if let Some(parent) = log_path.parent() {
        if !parent.as_os_str().is_empty() && parent.as_os_str() != std::ffi::OsStr::new("") {
            fs::create_dir_all(parent)
                .map_err(|e| anyhow::anyhow!("Cannot create log directory: {}", e))?;
        }
    }

    let log_file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(truncate)  // truncate=true to overwrite, truncate=false to append
        .append(!truncate)   // append=true if NOT truncating
        .open(log_file_path)
        .map_err(|e| anyhow::anyhow!("Cannot open log file: {}", e))?;

    env_logger::Builder::new()
        .format(|buf, record| {
            use chrono::Local;
            let timestamp = Local::now().format("%Y-%m-%dT%H:%M:%SZ");
            writeln!(
                buf,
                "[{} {:5} {}] {}",
                timestamp,
                record.level(),
                record.target(),
                record.args()
            )
        })
        .target(env_logger::Target::Pipe(Box::new(TeeWriter::new(log_file))))
        .filter_level(LevelFilter::Info)
        .try_init()
        .ok(); // Ignore if logger already initialized

    Ok(())
}

/// Normalize output path to end with /
pub fn normalize_path(path: &str) -> String {
    if path.ends_with('/') {
        path.to_string()
    } else {
        format!("{}/", path)
    }
}

/// Validate read datatype
pub fn validate_datatype(datatype: &str) -> Result<()> {
    match datatype {
        "clr" | "hifi" | "nanopore" => Ok(()),
        _ => anyhow::bail!("Invalid datatype. Must be one of: clr, hifi, nanopore. Got: {}", datatype),
    }
}

/// Get minimap2 preset based on read type
pub fn get_minimap2_preset(datatype: &str) -> &'static str {
    match datatype {
        "hifi" => "map-hifi",
        "nanopore" => "map-ont",
        "clr" | _ => "map-pb",
    }
}

/// Format bytes as human-readable size
pub fn format_size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut idx = 0;

    while size >= 1024.0 && idx < UNITS.len() - 1 {
        size /= 1024.0;
        idx += 1;
    }

    format!("{:.2} {}", size, UNITS[idx])
}

/// Count total number of reads (FASTQ records) in a gzipped FASTQ file.
/// Uses noodles for robust record parsing and MultiGzDecoder to handle
/// bgzipped (multi-block) files produced by sequencers.
pub fn count_fastq_reads(fastq_path: &str) -> Result<usize> {
    use flate2::read::MultiGzDecoder;
    use noodles::fastq::io::Reader as FastqReader;
    use std::io::BufReader;

    let file = fs::File::open(fastq_path)
        .map_err(|e| anyhow::anyhow!("Cannot open FASTQ file '{}': {}", fastq_path, e))?;
    let decoder = MultiGzDecoder::new(file);
    let mut reader = FastqReader::new(BufReader::new(decoder));

    let mut count = 0usize;
    for result in reader.records() {
        result.map_err(|e| anyhow::anyhow!("Error reading FASTQ record in '{}': {}", fastq_path, e))?;
        count += 1;
    }
    Ok(count)
}

#[derive(Debug, Clone, Copy, Default)]
pub struct FastqStats {
    pub reads: usize,
    pub bases: u64,
}

#[derive(Debug, Clone)]
pub struct SplitFastqResult {
    pub paths: Vec<String>,
    pub total_reads: usize,
    pub written_reads: usize,
    pub total_bases: u64,
    pub written_bases: u64,
    pub part_counts: Vec<usize>,
}

/// Count total reads and bases in a gzipped FASTQ file.
pub fn get_fastq_stats(fastq_path: &str) -> Result<FastqStats> {
    use flate2::read::MultiGzDecoder;
    use noodles::fastq::io::Reader as FastqReader;
    use std::io::BufReader;

    let file = fs::File::open(fastq_path)
        .map_err(|e| anyhow::anyhow!("Cannot open FASTQ file '{}': {}", fastq_path, e))?;
    let decoder = MultiGzDecoder::new(file);
    let mut reader = FastqReader::new(BufReader::new(decoder));

    let mut stats = FastqStats::default();
    for result in reader.records() {
        let record = result.map_err(|e| anyhow::anyhow!("Error reading FASTQ record in '{}': {}", fastq_path, e))?;
        stats.reads += 1;
        stats.bases += record.sequence().len() as u64;
    }

    Ok(stats)
}

/// Run seqkit stats on a FASTQ file and return read count and total bases.
/// Much faster than reading all records in Rust for large files.
pub fn seqkit_stats(fastq_path: &str) -> Result<FastqStats> {
    let output = std::process::Command::new("seqkit")
        .arg("stats")
        .arg("-j")
        .arg("8")  // use multiple threads for speed
        .arg("-T")   // tab-delimited output
        .arg(fastq_path)
        .output()
        .map_err(|e| anyhow::anyhow!("Failed to run seqkit stats on '{}': {}", fastq_path, e))?;

    if !output.status.success() {
        anyhow::bail!(
            "seqkit stats failed on '{}': {}",
            fastq_path,
            String::from_utf8_lossy(&output.stderr)
        );
    }

    // Tab-delimited columns: file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
    let stdout = String::from_utf8_lossy(&output.stdout);
    for line in stdout.lines().skip(1) {  // skip header
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 5 {
            let reads = fields[3].trim().parse::<usize>().unwrap_or(0);
            let bases = fields[4].trim().parse::<u64>().unwrap_or(0);
            return Ok(FastqStats { reads, bases });
        }
    }

    anyhow::bail!("Could not parse seqkit stats output for '{}'", fastq_path)
}

/// Split a FASTQ file into num_parts parts using seqkit split2.
/// Returns sorted list of output file paths.
pub fn seqkit_split2(input_path: &str, outdir: &str, num_parts: usize) -> Result<Vec<String>> {
    let status = std::process::Command::new("seqkit")
        .arg("split2")
        .arg("-p")
        .arg(num_parts.to_string())
        .arg("-O")
        .arg(outdir)
        .arg("-f")   // force overwrite existing files
        .arg(input_path)
        .status()
        .map_err(|e| anyhow::anyhow!("Failed to run seqkit split2 on '{}': {}", input_path, e))?;

    if !status.success() {
        anyhow::bail!("seqkit split2 failed on '{}'", input_path);
    }

    // seqkit names outputs <basename>.part_NNN.<ext> placed in outdir
    let mut paths: Vec<String> = fs::read_dir(outdir)
        .map_err(|e| anyhow::anyhow!("Cannot list split output dir '{}': {}", outdir, e))?
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| {
            let path = entry.path();
            let name = path.file_name()?.to_str()?.to_string();
            if name.contains(".part_") && (name.ends_with(".fastq.gz") || name.ends_with(".fq.gz")) {
                Some(path.to_string_lossy().to_string())
            } else {
                None
            }
        })
        .collect();

    paths.sort();

    if paths.is_empty() {
        anyhow::bail!(
            "seqkit split2 produced no output files in '{}' for input '{}'",
            outdir, input_path
        );
    }

    Ok(paths)
}

/// Subsample a FASTQ file to retain approximately `fraction` of reads using seqkit sample.
pub fn seqkit_sample(input_path: &str, output_path: &str, fraction: f64) -> Result<()> {
    let status = std::process::Command::new("seqkit")
        .arg("sample")
        .arg("-p")
        .arg(format!("{:.6}", fraction))
        .arg("-s")   // fixed seed for reproducibility
        .arg("42")
        .arg("-o")
        .arg(output_path)
        .arg(input_path)
        .status()
        .map_err(|e| anyhow::anyhow!("Failed to run seqkit sample on '{}': {}", input_path, e))?;

    if !status.success() {
        anyhow::bail!("seqkit sample failed on '{}'", input_path);
    }

    Ok(())
}

/// Split a gzipped FASTQ file into N parts and write to separate gzipped files.
/// Uses noodles for robust record parsing and MultiGzDecoder to handle bgzipped
/// (multi-block gzip) files produced by sequencers.
/// Returns paths plus input/output read/base counts.
pub fn split_fastq_gz(
    input_path: &str,
    outdir: &str,
    num_parts: usize,
    keep_fraction: Option<f64>,
) -> Result<SplitFastqResult> {
    use flate2::read::MultiGzDecoder;
    use flate2::write::GzEncoder;
    use noodles::fastq::io::{Reader as FastqReader, Writer as FastqWriter};
    use noodles::fastq::Record;
    use noodles::fastq::record::Definition;
    use std::io::{BufReader, BufWriter};

    let outdir = outdir.trim_end_matches('/');

    // Open input with MultiGzDecoder so both regular .gz and bgzipped files work
    let input_file = fs::File::open(input_path)
        .map_err(|e| anyhow::anyhow!("Cannot read FASTQ file '{}': {}", input_path, e))?;
    let decoder = MultiGzDecoder::new(input_file);
    let mut reader = FastqReader::new(BufReader::new(decoder));

    // Build output paths
    let paths: Vec<String> = (0..num_parts)
        .map(|i| format!("{}/split_{:02}.fastq.gz", outdir, i))
        .collect();

    // Create gzip writers — noodles fastq::Writer wraps any Write target
    let mut writers: Vec<FastqWriter<BufWriter<GzEncoder<fs::File>>>> = paths.iter()
        .map(|path| {
            let file = fs::File::create(path)
                .map_err(|e| anyhow::anyhow!("Cannot create split file '{}': {}", path, e))?;
            let encoder = GzEncoder::new(file, flate2::Compression::default());
            Ok(FastqWriter::new(BufWriter::new(encoder)))
        })
        .collect::<Result<Vec<_>>>()?;

    // Stream records and distribute round-robin
    let mut total = 0usize;
    let mut written_reads = 0usize;
    let mut total_bases = 0u64;
    let mut written_bases = 0u64;
    let mut current_part = 0usize;
    let mut part_counts = vec![0usize; num_parts];
    let keep_fraction = keep_fraction.unwrap_or(1.0).clamp(0.0, 1.0);

    for result in reader.records() {
        let record = result.map_err(|e| anyhow::anyhow!("Error reading FASTQ record: {}", e))?;
        total += 1;
        let record_bases = record.sequence().len() as u64;
        total_bases += record_bases;

        let keep_record = if keep_fraction >= 1.0 {
            true
        } else if keep_fraction <= 0.0 {
            false
        } else {
            let desired_written_bases = (total_bases as f64 * keep_fraction).floor() as u64;
            written_bases < desired_written_bases
        };

        if keep_record {
            writers[current_part].write_record(&record)
                .map_err(|e| anyhow::anyhow!("Error writing split FASTQ record: {}", e))?;
            part_counts[current_part] += 1;
            written_reads += 1;
            written_bases += record_bases;
            current_part = (current_part + 1) % num_parts;
        }
    }

    // Flush BufWriters and finish gzip encoders
    for writer in writers {
        let buf_writer = writer.into_inner();
        let gz_encoder = buf_writer.into_inner()
            .map_err(|e| anyhow::anyhow!("Error flushing split FASTQ writer: {}", e))?;
        gz_encoder.finish()
            .map_err(|e| anyhow::anyhow!("Error finalizing gzip split file: {}", e))?;
    }

    Ok(SplitFastqResult {
        paths,
        total_reads: total,
        written_reads,
        total_bases,
        written_bases,
        part_counts,
    })
}

/// Compute N50 statistic from a list of contig lengths
pub fn compute_n50(mut lengths: Vec<u64>) -> u64 {
    if lengths.is_empty() {
        return 0;
    }

    lengths.sort_by(|a, b| b.cmp(a)); // descending
    let total: u64 = lengths.iter().sum();
    let half_total = total / 2;

    let mut sum = 0;
    for len in lengths {
        sum += len;
        if sum >= half_total {
            return len;
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_path() {
        assert_eq!(normalize_path("./output"), "./output/");
        assert_eq!(normalize_path("./output/"), "./output/");
    }

    #[test]
    fn test_compute_n50() {
        let lengths = vec![100, 200, 300, 400, 500];
        let n50 = compute_n50(lengths);
        assert_eq!(n50, 400); // total=1500, half=750, sum reaches 750 at 400
    }

    #[test]
    fn test_validate_datatype() {
        assert!(validate_datatype("clr").is_ok());
        assert!(validate_datatype("hifi").is_ok());
        assert!(validate_datatype("nanopore").is_ok());
        assert!(validate_datatype("invalid").is_err());
    }

    #[test]
    fn test_get_minimap2_preset() {
        assert_eq!(get_minimap2_preset("hifi"), "map-hifi");
        assert_eq!(get_minimap2_preset("nanopore"), "map-ont");
        assert_eq!(get_minimap2_preset("clr"), "map-pb");
    }

    #[test]
    fn test_split_fastq_gz() {
        use flate2::write::GzEncoder;
        use noodles::fastq::io::Writer as FastqWriter;
        use noodles::fastq::Record;
        use noodles::fastq::record::Definition;
        use std::io::BufWriter;

        // Create a temporary directory
        let temp_dir = std::path::PathBuf::from("/tmp/inspector_test_split");
        let _ = std::fs::remove_dir_all(&temp_dir);
        std::fs::create_dir_all(&temp_dir).unwrap();

        // Create a test gzipped FASTQ file with 8 reads using noodles
        let input_fastq = temp_dir.join("test.fastq.gz");
        let names = ["read1","read2","read3","read4","read5","read6","read7","read8"];
        let seqs  = ["ACGT","TGCA","ACGTACGT","TTTT","AAAA","GGGG","CCCC","TAGC"];
        {
            let file = std::fs::File::create(&input_fastq).unwrap();
            let encoder = GzEncoder::new(file, flate2::Compression::default());
            let mut writer = FastqWriter::new(BufWriter::new(encoder));
            for (name, seq) in names.iter().zip(seqs.iter()) {
                let qual: Vec<u8> = vec![b'I'; seq.len()];
                let record = Record::new(
                    Definition::new(*name, ""),
                    seq.as_bytes().to_vec(),
                    qual,
                );
                writer.write_record(&record).unwrap();
            }
            let buf = writer.into_inner();
            let gz  = buf.into_inner().unwrap();
            gz.finish().unwrap();
        }

        // Split into 3 parts
        let output_dir = temp_dir.join("split_output");
        std::fs::create_dir_all(&output_dir).unwrap();

        let split_result = split_fastq_gz(
            input_fastq.to_str().unwrap(),
            output_dir.to_str().unwrap(),
            3,
            None,
        ).expect("split_fastq_gz failed");

        // Verify
        assert_eq!(split_result.total_reads, 8, "Expected 8 reads in total");
        assert_eq!(split_result.written_reads, 8, "Expected 8 written reads");
        assert_eq!(split_result.paths.len(), 3, "Expected 3 split files");

        // Verify each split file exists and check per-part counts returned directly (no re-read)
        assert_eq!(split_result.part_counts.len(), 3, "Expected 3 part count entries");
        // With 8 reads and 3 parts: 3, 3, 2 distribution (round-robin)
        assert_eq!(split_result.part_counts[0], 3, "Expected 3 reads in split part 0");
        assert_eq!(split_result.part_counts[1], 3, "Expected 3 reads in split part 1");
        assert_eq!(split_result.part_counts[2], 2, "Expected 2 reads in split part 2");

        for (i, file) in split_result.paths.iter().enumerate() {
            assert!(std::path::Path::new(file).exists(), "Split file {} does not exist", i);
        }

        // Cleanup
        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    #[test]
    fn test_count_fastq_reads() {
        use flate2::write::GzEncoder;
        use noodles::fastq::io::Writer as FastqWriter;
        use noodles::fastq::Record;
        use noodles::fastq::record::Definition;
        use std::io::BufWriter;

        // Create a temporary directory
        let temp_dir = std::path::PathBuf::from("/tmp/inspector_test_count");
        let _ = std::fs::remove_dir_all(&temp_dir);
        std::fs::create_dir_all(&temp_dir).unwrap();

        // Create a test gzipped FASTQ file with 5 reads using noodles
        let input_fastq = temp_dir.join("test.fastq.gz");
        {
            let file = std::fs::File::create(&input_fastq).unwrap();
            let encoder = GzEncoder::new(file, flate2::Compression::default());
            let mut writer = FastqWriter::new(BufWriter::new(encoder));
            for i in 1..=5 {
                let name = format!("read{}", i);
                let record = Record::new(
                    Definition::new(name.as_str(), ""),
                    b"ACGT".to_vec(),
                    b"IIII".to_vec(),
                );
                writer.write_record(&record).unwrap();
            }
            let buf = writer.into_inner();
            let gz  = buf.into_inner().unwrap();
            gz.finish().unwrap();
        }

        // Count reads
        let count = count_fastq_reads(input_fastq.to_str().unwrap())
            .expect("count_fastq_reads failed");

        assert_eq!(count, 5, "Expected 5 reads in test FASTQ file");

        let stats = get_fastq_stats(input_fastq.to_str().unwrap())
            .expect("get_fastq_stats failed");
        assert_eq!(stats.reads, 5, "Expected 5 reads in FASTQ stats");
        assert_eq!(stats.bases, 20, "Expected 20 bases in FASTQ stats");

        // Cleanup
        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    #[test]
    fn test_split_fastq_gz_with_subsampling() {
        use flate2::write::GzEncoder;
        use noodles::fastq::io::Writer as FastqWriter;
        use noodles::fastq::Record;
        use noodles::fastq::record::Definition;
        use std::io::BufWriter;

        let temp_dir = std::path::PathBuf::from("/tmp/inspector_test_split_subsample");
        let _ = std::fs::remove_dir_all(&temp_dir);
        std::fs::create_dir_all(&temp_dir).unwrap();

        let input_fastq = temp_dir.join("test.fastq.gz");
        {
            let file = std::fs::File::create(&input_fastq).unwrap();
            let encoder = GzEncoder::new(file, flate2::Compression::default());
            let mut writer = FastqWriter::new(BufWriter::new(encoder));
            for i in 1..=8 {
                let name = format!("read{}", i);
                let record = Record::new(
                    Definition::new(name.as_str(), ""),
                    b"ACGT".to_vec(),
                    b"IIII".to_vec(),
                );
                writer.write_record(&record).unwrap();
            }
            let buf = writer.into_inner();
            let gz = buf.into_inner().unwrap();
            gz.finish().unwrap();
        }

        let output_dir = temp_dir.join("split_output");
        std::fs::create_dir_all(&output_dir).unwrap();

        let split_result = split_fastq_gz(
            input_fastq.to_str().unwrap(),
            output_dir.to_str().unwrap(),
            4,
            Some(0.5),
        ).expect("split_fastq_gz failed");

        assert_eq!(split_result.total_reads, 8, "Expected 8 reads in total");
        assert_eq!(split_result.written_reads, 4, "Expected 4 reads after 0.5 subsampling");
        assert_eq!(split_result.total_bases, 32, "Expected 32 total bases");
        assert_eq!(split_result.written_bases, 16, "Expected 16 written bases after 0.5 subsampling");
        assert_eq!(split_result.part_counts.iter().sum::<usize>(), 4, "Expected 4 reads across split parts");

        // Cleanup
        let _ = std::fs::remove_dir_all(&temp_dir);
    }
}
