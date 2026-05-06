# Inspector-rs (Rust rewrite)

Inspector-rs is a rewrite of the reference-free long-read assembly evaluator and error-correction tool Inspector by Maggie Chen and others. 
Link: https://github.com/Maggi-Chen/Inspector (MIT license)

It provides:
- Assembly evaluation from long-read mappings
- Structural and small-scale error detection
- `long_read_QV` reporting
- Optional assembly correction workflow
- Logging to both screen and `Inspector.log`

## Motivation for rewrite
- Mapping can be very slow for large genomes - so this is sped up by dividing read sets into 4 parts
- Some python methods were unreliable or not robust for large read sets or genomes
- The original Inspector requires Python 2.7 and a very old minimap2, which are now way out of data and not allowed in our environment.

## Current Status

- Core evaluate/correct pipelines are implemented and runnable.
- Plotting functions are currently placeholders and not functional.

# Limitations
* This is not a 1:1 rewrite of Inspector, but should be very similar. We did not agree with some design choices in Inspector (eg allowing a base error to be predicted by just one read in a low coverage 5x dataset) so made some changes.
* Some parts may still be buggy - we are testing widely and welcome feedback.
* Plots are non-functional
* Collapses are found more rarely than the python implementation - may be buggy
* Evaluation should work, but the assembly correction mode is not of interest to me, is untested and may not work.

## Dependencies

## Download a release (make sure you get the most current one!)

```bash
wget https://github.com/colindaven/inspector-rs/releases/download/v0.0.4/inspector-rs
chmod a+x inspector-rs
# optional, or add to path in your .bashrc etc
sudo cp inspector-rs /usr/local/bin
cp inspector-rs ~/bin

```


## Rust - (only if modifying and compiling as a dev - else just download a release)

- Rust toolchain (edition 2021) 
- Cargo

Build:

```bash
cargo build --release
```

Binary path:

```bash
./target/release/inspector
```

## External tools (important)

Inspector shells out to these command-line tools - they must be in your PATH:

- `minimap2` (required): read-to-contig mapping
- `samtools` (required): BAM sort/merge/index/stats and mpileup
- `seqkit` (required): FASTQ stats, subsampling, splitting
- `flye` (optional): local re-assembly during correction for complex structural regions

Notes:
- The evaluate workflow performs preflight checks and logs tool versions to screen and `Inspector.log`.
- If `flye` is unavailable, correction falls back to simple structural patching logic.

## What Inspector Does

## 1 Evaluate assembly quality - part implemented

High-level steps:
- Validate and filter contigs
- Compute assembly stats (N50, lengths)
- Estimate read coverage from FASTQ base counts
- Optional read subsampling (`--read-coverage`, default 30x)
- Read mapping to contigs
- Structural error signal detection and clustering
- Base-level error detection from pileup
- Compute and report `long_read_QV`

Outputs include:
- `valid_contig.fa`
- `summary_statistics`
- `read_to_contig.bam` (+ index)
- error workspaces (`ae_merge_workspace`, `base_error_workspace`, etc.)
- `Inspector.log`

## 2 Correct assembly errors - NOT IMPLEMENTED YET 

- Reads error predictions from evaluate output
- Applies base and structural corrections
- Optionally invokes Flye for local re-assembly attempts
- Writes corrected assembly and correction summary

## FASTQ handling and performance behavior

- Coverage is estimated as:

  `base_coverage = total_input_read_bp / genome_size_bp`

- Read subsampling behavior:
  - If estimated coverage > `--read-coverage`, reads are subsampled.
  - If estimated coverage cannot be calculated, no subsampling is performed.

- Split policy based on genome size:
  - If genome size > 500 Mbp: split each FASTQ into 4 parts.
  - If genome size <= 500 Mbp: do not split FASTQ files.

## Usage

Show help:

```bash
./target/release/inspector --help
./target/release/inspector evaluate --help
./target/release/inspector correct --help
./target/release/inspector all --help
```

## Evaluate

```bash
./target/release/inspector evaluate \
  --contig assembly.fasta \
  --read reads.fastq.gz \
  --datatype hifi \
  --outpath ./inspector-out/ \
  --thread 16
```

With explicit coverage target and depth:

```bash
./target/release/inspector evaluate \
  --contig assembly.fasta \
  --read reads.fastq.gz \
  --datatype nanopore_1041 \
  --outpath ./inspector-out/ \
  --thread 24 \
  --read-coverage 30 \
  --min-depth 10
```

Skip selected steps:

```bash
./target/release/inspector evaluate \
  --contig assembly.fasta \
  --read reads.fastq.gz \
  --datatype nanopore_1041 \
  --outpath ./inspector-out/ \
  --skip-read-mapping \
  --skip-structural-error \
  --skip-base-error
```

## Correct

```bash
./target/release/inspector correct \
  --input ./inspector-out/ \
  --assembly ./inspector-out/valid_contig.fa \
  --read reads.fastq.gz \
  --datatype hifi \
  --outpath ./inspector-correct-out/ \
  --thread 16 \
  --min-correction-support 1
```

## All (evaluate + correct)

```bash
./target/release/inspector all \
  --contig assembly.fasta \
  --read reads.fastq.gz \
  --datatype hifi \
  --outpath ./inspector-run/ \
  --thread 16 \
  --read-coverage 30 \
  --min-correction-support 1
```

---

## Parameters

### `evaluate` subcommand

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--contig` | `-c` | FILE (1+) | *(required)* | Assembly contigs in FASTA format. |
| `--read` | `-r` | FILE (1+) | *(required)* | Sequencing reads in FASTA/FASTQ format (can be gzipped). |
| `--datatype` | `-d` | TYPE | `nanopore_1041` | Input read type. Accepted values: `clr`, `hifi`, `nanopore_94`, `nanopore_1041`. |
| `--outpath` | `-o` | PATH | `./inspector-out/` | Output directory. |
| `--reference` | | FILE | *(none)* | Optional reference genome in FASTA format. |
| `--thread` | `-t` | INT | `8` | Number of threads. |
| `--read-coverage` | | INT | `30` | Target read coverage for pre-mapping subsampling. Reads are downsampled if estimated coverage exceeds this value. |
| `--min-depth` | | INT | *(auto: 20% of avg depth)* | Minimal read-alignment depth for a base to be included in QV calculation. |
| `--min-contig-length` | | INT | `10000` | Minimal contig length to be included in evaluation. |
| `--min-contig-length-assemblyerror` | | INT | `1000000` | Minimal contig length for structural error detection. |
| `--min-assembly-error-size` | | INT | `50` | Minimal size (bp) of assembly errors to report. |
| `--max-assembly-error-size` | | INT | `4000000` | Maximal size (bp) of assembly errors to report. |
| `--noplot` | | flag | off | Disable plot generation. |
| `--skip-read-mapping` | | flag | off | Skip the read-to-contig mapping step (use pre-existing BAM). |
| `--skip-structural-error` | | flag | off | Skip structural error clustering and filtering. |
| `--skip-structural-error-detect` | | flag | off | Skip structural error signal detection (also skips mapping depth collection). |
| `--skip-base-error` | | flag | off | Skip small-scale error counting and QV contribution. |
| `--skip-base-error-detect` | | flag | off | Skip pileup-based base error detection (use pre-existing pileup files). |

### `correct` subcommand

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--input` | `-i` | PATH | *(required)* | Inspector `evaluate` output directory containing error BED files. |
| `--read` | `-r` | FILE (1+) | *(required)* | Sequencing reads in FASTA/FASTQ format. |
| `--assembly` | `-a` | FILE | *(required)* | Assembly contigs in FASTA format (typically `valid_contig.fa` from evaluate output). |
| `--outpath` | `-o` | PATH | `./inspector-correct-out/` | Output directory for corrected assembly. |
| `--datatype` | `-d` | TYPE | `nanopore_1041` | Input read type. Accepted values: `clr`, `hifi`, `nanopore_94`, `nanopore_1041`. |
| `--thread` | `-t` | INT | `8` | Number of threads. |
| `--flye-timeout` | | INT | `3600` | Timeout in seconds for each Flye local-reassembly call. |
| `--min-correction-support` | | INT | `1` | Minimum supporting read count in BED file required to apply a correction. |

### `all` subcommand

Runs `evaluate` followed by `correct` in sequence. Evaluate output goes to `<outpath>/evaluate/`; correction output goes to `<outpath>/correct/`.

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--contig` | `-c` | FILE (1+) | *(required)* | Assembly contigs in FASTA format. |
| `--read` | `-r` | FILE (1+) | *(required)* | Sequencing reads in FASTA/FASTQ format (can be gzipped). |
| `--datatype` | `-d` | TYPE | `nanopore_1041` | Input read type. Accepted values: `clr`, `hifi`, `nanopore`. |
| `--outpath` | `-o` | PATH | `./inspector-out/` | Root output directory. |
| `--reference` | | FILE | *(none)* | Optional reference genome in FASTA format. |
| `--thread` | `-t` | INT | `8` | Number of threads. |
| `--read-coverage` | | INT | `30` | Target read coverage for pre-mapping subsampling. |
| `--min-depth` | | INT | *(auto: 20% of avg depth)* | Minimal read-alignment depth for QV calculation. |
| `--min-contig-length` | | INT | `10000` | Minimal contig length to be included in evaluation. |
| `--min-contig-length-assemblyerror` | | INT | `1000000` | Minimal contig length for structural error detection. |
| `--min-assembly-error-size` | | INT | `50` | Minimal size (bp) of assembly errors to report. |
| `--max-assembly-error-size` | | INT | `4000000` | Maximal size (bp) of assembly errors to report. |
| `--noplot` | | flag | off | Disable plot generation. |
| `--skip-read-mapping` | | flag | off | Skip the read-to-contig mapping step. |
| `--skip-structural-error` | | flag | off | Skip structural error detection and clustering. |
| `--skip-base-error` | | flag | off | Skip small-scale error detection. |
| `--flye-timeout` | | INT | `3600` | Timeout in seconds for each Flye local-reassembly call. |
| `--min-correction-support` | | INT | `1` | Minimum supporting read count required to apply a correction. |

---

## Logging

- Logs are written to both terminal and `Inspector.log` in the chosen output directory.
- Log file is overwritten at start of each run.
- Preflight tool checks and versions are logged.

## Datatype values

Supported `--datatype` values:
- `nanopore_1041` - Oxford Nanopore reads (chemistry 10.4.1 or newer, higher accuracy) → uses `lr:hq` preset
- `nanopore_94` - Oxford Nanopore reads (chemistry 9.4 or older, noisier reads) → uses `map-ont` preset
- `hifi` - PacBio HiFi (High Fidelity) reads
- `clr` - PacBio CLR (Continuous Long Reads)

## Clustering algorithm

The clustering algorithm processes structural variant signals to identify and reconcile assembly errors through a sophisticated multi-stage pipeline involving signal grouping, genotyping, and error filtering.

### Key Parameters and Default Values

#### Coverage-Based Thresholds
- **`coverage`**: Average read coverage calculated from BAM mapping statistics (default: 20 if calculation fails)
- **`highcov`**: `coverage × 2` (default: 40 when coverage=20)
- **`lowcov`**: `coverage ÷ 2` (default: 10 when coverage=20)

#### Additional Parameters
- **`min_support`**: Minimum supporting reads (typically 3-5)
- **`max_depth`**: Maximum depth for small SV clustering (`coverage × 2`)
- **`rat`**: Coverage ratio threshold (0.8 for HiFi, 0.7 for other datatypes)
- **`min_size`**: Minimum SV size threshold for reporting

### Algorithm Stages

#### Stage 1: Signal Clustering

The algorithm processes three types of structural variants separately:

##### 1.1 Deletion Clustering (`cluster`)
- **Large deletions** (>2000 bp): Adaptive window clustering with 1600bp windows
- **Small deletions** (≤3000 bp): Depth-map spatial clustering
  - Creates depth array across contig length
  - Marks deletion spans: `[pos-1, pos+size)`
  - Identifies high-depth regions (>3× threshold) with max depth ≤ `max_depth`
  - Groups signals within these regions

##### 1.2 Insertion/Inversion Clustering (`cluster_insertions`)
- **Large insertions/inversions** (>2000 bp): Adaptive window clustering
- **Small insertions/inversions** (≤3000 bp): Position-based clustering
  - Marks 202bp windows: `[pos-101, pos+101]`
  - Assigns signals within 100bp of region boundaries
  - Handles point-like events (insertions don't extend genomically)

##### 1.3 Deduplication
- Removes large SVs that significantly overlap with small SVs
- Prevents double-counting of the same structural variant

#### Stage 2: Genotyping (`genotype`)

For each clustered SV:
1. **Flanking coverage calculation**: Count reads in 100bp windows before and after SV
2. **Genotype assignment**:
   - `1/1` (homozygous): `supporting_reads ≥ 0.6 × min(left_cov, right_cov)`
   - `1/0` (heterozygous): Otherwise

#### Stage 3: Error Filtering and Reconciliation (`filter_errors`)

This is where **`highcov`** and **`lowcov`** play their crucial roles:

##### 3.1 Expansion/Collapse Reconciliation
The algorithm attempts to match expansion and collapse signals that represent the same underlying assembly error:

- **Spatial proximity**: Collapse start must be within [expansion_start-250, expansion_end+250]
- **Size compatibility**: Expansion size must be < 20× collapse size
- **Coverage ratio analysis**: `expansion_support / collapse_support` must be in [0.33, 3.0]

**Reconciliation outcomes**:
- **Balanced evidence** (ratio 0.33-3.0): Create specific Expansion or Collapse events
- **Expansion dominant** (ratio >3.0): Report as Expansion 
- **Collapse dominant** (ratio <0.33): Report as Collapse
- **Haplotype switch**: When total unique reads exceed individual thresholds

##### 3.2 Final Quality Filtering

Each SV must pass **all** of these criteria:

1. **Minimum support**: `supporting_reads ≥ 10`
2. **Coverage ratio**: `supporting_reads ≥ rat × depth_min`
   - HiFi: `supporting_reads ≥ 0.8 × depth_min`  
   - Other: `supporting_reads ≥ 0.7 × depth_min`
3. **Coverage bounds**: `lowcov ≤ depth_min < highcov`
4. **Size threshold**: `max_size ≥ min_size`

### How `highcov` and `lowcov` Affect Results

#### Expansion Detection
- **`lowcov` threshold**: Ensures sufficient baseline coverage to distinguish real expansions from low-coverage artifacts
  - Too low → false positives in low-coverage regions
  - Too high → missed expansions in moderately covered regions

- **`highcov` threshold**: Prevents calling expansions in extremely high-coverage regions where read pileup might create false signals
  - Too low → missed expansions in high-coverage regions
  - Too high → false positives from repetitive regions

#### Collapse Detection  
- **Coverage bounds work similarly**: Collapses need sufficient flanking coverage to be confidently called
- **Depth filtering**: Prevents calling collapses in regions with aberrant coverage patterns

#### Practical Impact

**With default coverage=20**:
- `lowcov=10, highcov=40`
- **Optimal range**: SVs with flanking depth 10-39× are candidates
- **Filtered out**: SVs in very low coverage (<10×) or very high coverage (≥40×) regions

**Effect of changing coverage**:
- **Higher coverage input** → higher thresholds → more stringent filtering → fewer but higher-confidence calls
- **Lower coverage input** → lower thresholds → more permissive filtering → more calls but potentially more false positives

The algorithm balances sensitivity and specificity by using these coverage-based thresholds to ensure structural variant calls are made only in regions with appropriate and stable read depth characteristics.


# Credits

Inspector-rs is a rewrite of the reference-free long-read assembly evaluator and error-correction tool Inspector by Maggie Chen and others - so thanks to the original authors for a nice approach!

Link: https://github.com/Maggi-Chen/Inspector (MIT license)
