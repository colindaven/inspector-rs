/// SV merging and clustering operations
/// Position-window based clustering of SV signals

pub mod merge_ops;
pub mod clustering;

pub use merge_ops::{merge_insertions, merge_deletions, merge_translocations};
pub use clustering::{cluster, cluster_insertions, genotype, filter_errors};

use anyhow::Result;

/// Main merge function coordinator
pub fn merge_all(
    outpath: &str,
    min_support: usize,
    datatype: &str,
) -> Result<()> {
    // Placeholder
    Ok(())
}
