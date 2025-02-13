#![doc = include_str!("../README.md")]

pub mod builds;
mod genome;

pub use genome::{Contig, GenomeBuild, GenomeBuildIdentifier};
