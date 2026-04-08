# Inspector (Rust rewrite)

Inspector is a reference-free long-read assembly evaluator and error-correction tool.

It provides:
- Assembly evaluation from long-read mappings
- Structural and small-scale error detection
- `long_read_QV` reporting
- Optional assembly correction workflow
- Logging to both screen and `Inspector.log`

## Motivation for rewrite
- Mapping can be very slow for large genomes - so this is sped up by dividing read sets into 4 parts
- Some python methods were unreliable or not robust for large read sets or genomes

## Current Status

- Core evaluate/correct pipelines are implemented and runnable.
- Plotting functions are currently placeholders and not functional.
  - `plot_n50`
  - `plot_dotplot`
  - `plot_na50`

## Dependencies

## Download a release

```bash
wget https://github.com/colindaven/inspector-rs/releases/download/v0.0.1/inspector
chmod a+x inspector
# optional, or add to path
sudo cp inspector /usr/local/bin
cp inspector ~/bin

```


## Rust - (only if compiling - else just download a release)

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

Inspector shells out to these command-line tools:

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
  --datatype nanopore \
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
  --datatype clr \
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

## Logging

- Logs are written to both terminal and `Inspector.log` in the chosen output directory.
- Log file is overwritten at start of each run.
- Preflight tool checks and versions are logged.

## Datatype values

Supported `--datatype` values:
- `clr`
- `hifi`
- `nanopore`

## Limitations / Notes

- Plot outputs are currently non-functional placeholders.
- Some modules still contain placeholder logic intended for further refinement.
- Performance and behavior depend strongly on external tool versions and input data characteristics.
