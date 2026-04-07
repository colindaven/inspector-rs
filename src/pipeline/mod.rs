/// Main pipeline orchestrators for evaluate and correct subcommands

pub mod evaluate;
pub mod correct;

pub use evaluate::{evaluate, EvaluateConfig};
pub use correct::{correct, CorrectConfig};
