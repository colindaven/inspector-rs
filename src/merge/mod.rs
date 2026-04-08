/// SV merging and clustering operations
/// Position-window based clustering of SV signals

pub mod merge_ops;
pub mod clustering;

pub use merge_ops::{merge_insertions, merge_deletions, merge_translocations};
pub use merge_ops::{merge_one_event, merge_with_bimodal};
pub use clustering::{cluster, cluster_insertions, genotype, filter_errors, write_structural_error_tsv, write_summary_statistics_extended};

use anyhow::Result;

/// Main merge function coordinator
pub fn merge_all(
    _outpath: &str,
    _min_support: usize,
    _datatype: &str,
) -> Result<()> {
    // Placeholder
    Ok(())
}
