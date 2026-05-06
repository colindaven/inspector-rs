use clap::{Parser, Subcommand};
use anyhow::Result;
use log::info;

#[derive(Parser)]
#[command(name = "inspector")]
#[command(about = "Reference-free assembly evaluator and error corrector", long_about = None, version = "1.0.0")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Evaluate assembly quality using long reads
    Evaluate {
        /// Assembly contigs in FASTA format (can specify 1 or 2 files for haploid/diploid)
        #[arg(short, long, value_name = "FILE", num_args = 1.., required = true)]
        contig: Vec<String>,

        /// Sequencing reads in FASTA/FASTQ format (can be gzipped)
        #[arg(short, long, value_name = "FILE", num_args = 1.., required = true)]
        read: Vec<String>,

        /// Input read type
        #[arg(short, long, value_name = "TYPE", default_value = "nanopore_1041")]
        datatype: String,

        /// Output directory
        #[arg(short, long, value_name = "PATH", default_value = "./inspector-out/")]
        outpath: String,

        /// Optional reference genome in FASTA format
        #[arg(long, value_name = "FILE")]
        reference: Option<String>,

        /// Number of threads
        #[arg(short, long, default_value = "8")]
        thread: usize,

        /// Target read coverage for pre-mapping subsampling; reads above this estimated coverage are downsampled.
        /// Set to 0 to disable subsampling (default: disabled for performance)
        #[arg(long, default_value = "0")]
        read_coverage: usize,

        /// Minimal read-alignment depth for a base to be considered in QV calculation.
        /// If not specified, defaults to 20% of average depth
        #[arg(long, value_name = "INT")]
        min_depth: Option<usize>,

        /// Minimal read mapping quality used for small-scale base-error detection
        #[arg(long, value_name = "INT", default_value = "20")]
        min_mapping_quality: u8,

        /// Minimal length for a contig to be evaluated
        #[arg(long, default_value = "10000")]
        min_contig_length: usize,

        /// Minimal contig length for assembly error detection
        #[arg(long, default_value = "1000000")]
        min_contig_length_assemblyerror: usize,

        /// Minimal size for assembly errors
        #[arg(long, default_value = "50")]
        min_assembly_error_size: usize,

        /// Maximal size for assembly errors
        #[arg(long, default_value = "4000000")]
        max_assembly_error_size: usize,

        /// Do not make plots
        #[arg(long)]
        noplot: bool,

        /// Skip the step of mapping reads to contig
        #[arg(long)]
        skip_read_mapping: bool,

        /// Skip the step of identifying large structural errors
        #[arg(long)]
        skip_structural_error: bool,

        /// Skip the step of detecting large structural errors
        #[arg(long)]
        skip_structural_error_detect: bool,

        /// Skip the step of identifying small-scale errors
        #[arg(long)]
        skip_base_error: bool,

        /// Skip the step of detecting small-scale errors from pileup
        #[arg(long)]
        skip_base_error_detect: bool,
    },

    /// Correct assembly errors using structural/base error predictions
    Correct {
        /// Inspector evaluation output directory
        #[arg(short, long, value_name = "PATH", required = true)]
        input: String,

        /// Sequencing reads in FASTA/FASTQ format
        #[arg(short, long, value_name = "FILE", num_args = 1.., required = true)]
        read: Vec<String>,

        /// Assembly contigs in FASTA format
        #[arg(short, long, value_name = "FILE", required = true)]
        assembly: String,

        /// Output directory for corrected assembly
        #[arg(short, long, value_name = "PATH", default_value = "./inspector-correct-out/")]
        outpath: String,

        /// Input read type
        #[arg(short, long, value_name = "TYPE", default_value = "nanopore_1041")]
        datatype: String,

        /// Number of threads
        #[arg(short, long, default_value = "8")]
        thread: usize,

        /// Timeout for Flye assembly (seconds)
        #[arg(long, default_value = "3600")]
        flye_timeout: u64,

        /// Minimal BED support (read count) required to apply an error correction
        #[arg(long, default_value = "1")]
        min_correction_support: usize,
    },

    /// Run full pipeline: evaluate assembly then correct errors (default mode)
    All {
        /// Assembly contigs in FASTA format (can specify 1 or 2 files for haploid/diploid)
        #[arg(short, long, value_name = "FILE", num_args = 1.., required = true)]
        contig: Vec<String>,

        /// Sequencing reads in FASTA/FASTQ format (can be gzipped)
        #[arg(short, long, value_name = "FILE", num_args = 1.., required = true)]
        read: Vec<String>,

        /// Input read type
        #[arg(short, long, value_name = "TYPE", default_value = "nanopore_1041")]
        datatype: String,

        /// Output directory (evaluation output will be in <outpath>/evaluate/,
        /// corrected assembly in <outpath>/correct/)
        #[arg(short, long, value_name = "PATH", default_value = "./inspector-out/")]
        outpath: String,

        /// Optional reference genome in FASTA format
        #[arg(long, value_name = "FILE")]
        reference: Option<String>,

        /// Number of threads
        #[arg(short, long, default_value = "8")]
        thread: usize,

        /// Target read coverage for pre-mapping subsampling; reads above this estimated coverage are downsampled.
        /// Set to 0 to disable subsampling (default: disabled for performance)
        #[arg(long, default_value = "0")]
        read_coverage: usize,

        /// Minimal read-alignment depth for a base to be considered in QV calculation
        #[arg(long, value_name = "INT")]
        min_depth: Option<usize>,

        /// Minimal read mapping quality used for small-scale base-error detection
        #[arg(long, value_name = "INT", default_value = "20")]
        min_mapping_quality: u8,

        /// Minimal length for a contig to be evaluated
        #[arg(long, default_value = "10000")]
        min_contig_length: usize,

        /// Minimal contig length for assembly error detection
        #[arg(long, default_value = "1000000")]
        min_contig_length_assemblyerror: usize,

        /// Minimal size for assembly errors
        #[arg(long, default_value = "50")]
        min_assembly_error_size: usize,

        /// Maximal size for assembly errors
        #[arg(long, default_value = "4000000")]
        max_assembly_error_size: usize,

        /// Do not make plots
        #[arg(long)]
        noplot: bool,

        /// Skip the step of mapping reads to contig
        #[arg(long)]
        skip_read_mapping: bool,

        /// Skip the step of identifying large structural errors
        #[arg(long)]
        skip_structural_error: bool,

        /// Skip the step of identifying small-scale errors
        #[arg(long)]
        skip_base_error: bool,

        /// Timeout for Flye assembly correction (seconds)
        #[arg(long, default_value = "3600")]
        flye_timeout: u64,

        /// Minimal BED support (read count) required to apply an error correction
        #[arg(long, default_value = "1")]
        min_correction_support: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Evaluate {
            contig,
            read,
            datatype,
            outpath,
            reference,
            thread,
            read_coverage,
            min_depth,
            min_mapping_quality,
            min_contig_length,
            min_contig_length_assemblyerror,
            min_assembly_error_size,
            max_assembly_error_size,
            noplot,
            skip_read_mapping,
            skip_structural_error,
            skip_structural_error_detect,
            skip_base_error,
            skip_base_error_detect,
        } => {
            // Initialize file logger FIRST (before any log messages)
            // Log file goes in the output root directory (not subdir), truncate on start
            let normalized_outpath = inspector::utils::normalize_path(&outpath);
            let log_file_path = format!("{}Inspector.log", normalized_outpath);
            inspector::utils::init_file_logger(&log_file_path, true)
                .unwrap_or_else(|e| eprintln!("Warning: Could not setup file logging: {}", e));

            info!("Starting assembly evaluation");
            info!("Contigs: {:?}", contig);
            info!("Reads: {:?}", read);
            info!("Datatype: {}", datatype);
            info!("Output: {}", outpath);
            info!("Threads: {}", thread);

            inspector::evaluate::run(inspector::evaluate::EvaluateConfig {
                contig,
                read,
                datatype,
                outpath,
                reference,
                thread,
                read_coverage,
                min_depth,
                min_mapping_quality,
                min_contig_length,
                min_contig_length_assemblyerror,
                min_assembly_error_size,
                max_assembly_error_size,
                noplot,
                skip_read_mapping,
                skip_structural_error,
                skip_structural_error_detect,
                skip_base_error,
                skip_base_error_detect,
            })?;
        }

        Commands::Correct {
            input,
            read,
            assembly,
            outpath,
            datatype,
            thread,
            flye_timeout,
            min_correction_support,
        } => {
            // Initialize file logger for correct command, truncate on start
            let normalized_outpath = inspector::utils::normalize_path(&outpath);
            let log_file_path = format!("{}Inspector.log", normalized_outpath);
            inspector::utils::init_file_logger(&log_file_path, true)
                .unwrap_or_else(|e| eprintln!("Warning: Could not setup file logging: {}", e));
            info!("Starting assembly correction");
            info!("Input: {}", input);
            info!("Assembly: {}", assembly);
            info!("Output: {}", outpath);

            inspector::correct::run(inspector::correct::CorrectConfig {
                input,
                read,
                assembly,
                outpath,
                datatype,
                thread,
                flye_timeout,
                min_correction_support,
            })?;
        }

        Commands::All {
            contig,
            read,
            datatype,
            outpath,
            reference,
            thread,
            read_coverage,
            min_depth,
            min_mapping_quality,
            min_contig_length,
            min_contig_length_assemblyerror,
            min_assembly_error_size,
            max_assembly_error_size,
            noplot,
            skip_read_mapping,
            skip_structural_error,
            skip_base_error,
            flye_timeout,
            min_correction_support,
        } => {
            use inspector::utils::normalize_path;
            let outpath = normalize_path(&outpath);
            
            // Initialize file logger for all command at root level, truncate on start
            let log_file_path = format!("{}Inspector.log", outpath);
            inspector::utils::init_file_logger(&log_file_path, true)
                .unwrap_or_else(|e| eprintln!("Warning: Could not setup file logging: {}", e));
            
            let eval_outpath = format!("{}evaluate/", outpath);
            let correct_outpath = format!("{}correct/", outpath);

            // Phase 1: Evaluate
            info!("=== Inspector: Step 1/2 — Assembly Evaluation ===");
            info!("Contigs: {:?}", contig);
            info!("Reads: {:?}", read);
            info!("Datatype: {}", datatype);
            info!("Evaluation output: {}", eval_outpath);
            info!("Threads: {}", thread);

            inspector::evaluate::run(inspector::evaluate::EvaluateConfig {
                contig: contig.clone(),
                read: read.clone(),
                datatype: datatype.clone(),
                outpath: eval_outpath.clone(),
                reference,
                thread,
                read_coverage,
                min_depth,
                min_mapping_quality,
                min_contig_length,
                min_contig_length_assemblyerror,
                min_assembly_error_size,
                max_assembly_error_size,
                noplot,
                skip_read_mapping,
                skip_structural_error,
                skip_structural_error_detect: false,
                skip_base_error,
                skip_base_error_detect: false,
            })?;

            // The assembly file used is valid_contig.fa produced by evaluation
            let corrected_assembly_input = format!("{}valid_contig.fa", eval_outpath);

            // Phase 2: Correct
            info!("=== Inspector: Step 2/2 — Assembly Correction ===");
            info!("Reading errors from: {}", eval_outpath);
            info!("Correction output: {}", correct_outpath);

            inspector::correct::run(inspector::correct::CorrectConfig {
                input: eval_outpath,
                read,
                assembly: corrected_assembly_input,
                outpath: correct_outpath,
                datatype,
                thread,
                flye_timeout,
                min_correction_support,
            })?;

            info!("=== Inspector complete ===");
        }
    }

    Ok(())
}
