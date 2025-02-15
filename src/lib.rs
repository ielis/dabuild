//! # dabuild
//!
//! `dabuild` provides you with genome build metadata.
//!
//! ## Examples
//!
//! ### Load genome build
//!
//! The [`builds`] module provides several bundled builds.
//! Alternatively, you can load a build from an assembly report.
//!
//! See the [`builds`] documentation for more info.
//!
//! ### Use genome build
//!
//! Genome build is basically a data container and the usage involves accessing the data.
//!
//! ## Examples
//!
//! We show several examples with the *GRCh38.p13* genome build.
//!
//! ```rust
//! use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! use dabuild::builds::get_grch38_p13;
//!
//! let build: GenomeBuild<u32> = get_grch38_p13();
//! ```
//!
//! ### Check build identifiers
//!
//! We can check the major assembly and the patch of the build:
//!
//! ```rust
//! # use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! # use dabuild::builds::get_grch38_p13;
//! # let build: GenomeBuild<u32> = get_grch38_p13();
//!
//! assert_eq!(build.id().major_assembly(), "GRCh38");
//! assert_eq!(build.id().patch(), Some("p13"));
//! ```
//!
//! ### Access contigs
//!
//! The genome build contains one or more contigs.
//!
//! We can iterate over all contigs, e.g. to count them:
//!
//! ```rust
//! # use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! # use dabuild::builds::get_grch38_p13;
//! # let build: GenomeBuild<u32> = get_grch38_p13();
//!
//! let count = build.contigs().count();
//! assert_eq!(count, 640);
//! ```
//!
//! and we can also access a contig (e.g. `chrY`) by one of its names:
//!
//! ```rust
//! # use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! # use dabuild::builds::get_grch38_p13;
//! # let build: GenomeBuild<u32> = get_grch38_p13();
//!
//! // Query by name ...
//! let y = build.contig_by_name("Y");
//! assert!(y.is_some());
//!
//! /// ... or GenBank accession ...
//! let y = build.contig_by_name("CM000686.2");
//! assert!(y.is_some());
//!
//! /// ... or RefSeq accession ...
//! let y = build.contig_by_name("NC_000024.10");
//! assert!(y.is_some());
//!
//! /// ... or UCSC identifier.
//! let y = build.contig_by_name("chrY");
//! assert!(y.is_some());
//! ```

pub mod builds;
mod genome;

pub use genome::{Contig, GenomeBuild, GenomeBuildIdentifier};
