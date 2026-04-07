pub mod detect;
pub mod merge;
pub mod base_error;
pub mod correction;
pub mod plot;
pub mod static_analysis;
pub mod pipeline;
pub mod models;
pub mod utils;

pub use pipeline::{evaluate, correct};
