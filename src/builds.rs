//! Get the most commonly used genome builds.
//!
//! We provide several common genome build out-of-shelf.
//! Other builds can be loaded from a genome assembly report.
//!
//! ## Bundled genome builds
//!
//! These bundled genome builds can be loaded using the respective loader function:
//! * *GRCh37.p13*: [`get_grch37_p13`]
//! * *GRCh38.p13*: [`get_grch38_p13`]
//!
//! ### Example
//!
//! Load *GRCh38.p13* (*Homo sapiens*):
//!
//! ```rust
//! use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! use dabuild::builds::get_grch38_p13;
//!
//! let build: GenomeBuild = get_grch38_p13();
//! ```
//!
//! ## Load from an assembly report
//!
//! A genome build can be loaded from the Genome Reference Consortium assembly report
//! using [`parse_assembly_report`] function.
//!
//! For instance, *GRCm39* can be loaded from
//! an [example assembly report](https://github.com/ielis/dabuild/blob/master/data/GCF_000001635.27_GRCm39_assembly_report.txt):
//!
//! ```rust
//! # use dabuild::{GenomeBuild, GenomeBuildIdentifier};
//! # use dabuild::builds::parse_assembly_report;
//! use std::{fs::File, io::BufReader, str::FromStr};
//!
//! let path = "data/GCF_000001635.27_GRCm39_assembly_report.txt";
//! let build: GenomeBuild = parse_assembly_report(
//!         GenomeBuildIdentifier::from_str("GRCm39").expect("Infallible"),
//!         BufReader::new(File::open(path).expect("File should be present and readable")),
//! ).expect("No I/O or format issues");
//!
//! assert_eq!(build.id().major_assembly(), "GRCm39");
//! ```
//!

use super::{Contig, GenomeBuild, GenomeBuildIdentifier};
use std::num::ParseIntError;
use std::{io, io::BufRead};
use thiserror::Error;

#[allow(non_upper_case_globals)]
const GRCh37_p13: &[u8] = include_bytes!("data/GCF_000001405.25_GRCh37.p13_assembly_report.tsv");
#[allow(non_upper_case_globals)]
const GRCh38_p13: &[u8] = include_bytes!("data/GCF_000001405.39_GRCh38.p13_assembly_report.tsv");

/// Get the *GRCh37.p13* build.
///
/// ## Panics
///
/// If the builtin assembly report cannot be parsed (should not happen).
pub fn get_grch37_p13() -> GenomeBuild {
    let id = GenomeBuildIdentifier::from(("GRCh37", "p13"));
    parse_assembly_report(id, GRCh37_p13)
        .expect("The embedded assembly report for GRCh37.p13 should be valid")
}

/// Get the *GRCh38.p13* build.
///
/// ## Panics
///
/// If the builtin assembly report cannot be parsed (should not really happen).
pub fn get_grch38_p13() -> GenomeBuild {
    let id = GenomeBuildIdentifier::from(("GRCh38", "p13"));
    parse_assembly_report(id, GRCh38_p13)
        .expect("The embedded assembly report for GRCh38.p13 should be valid")
}

/// Parse an assembly report into a [`GenomeBuild`].
///
/// The assembly report is expected to include a header lines that start with `#`
/// and a tab-separated lines, one contig per line.
/// Each contig line is expected to contain the following 10 fields:
///
/// * Sequence-Name
/// * Sequence-Role
/// * Assigned-Molecule
/// * Assigned-Molecule-Location/Type
/// * GenBank-Accn
/// * Relationship
/// * RefSeq-Accn
/// * Assembly-Unit
/// * Sequence-Length
/// * UCSC-style-name
///
/// ## Errors
///
/// The parsing can fail from the reasons outlined in [`GenomeAssemblyParseError`].
///
pub fn parse_assembly_report<R>(
    id: GenomeBuildIdentifier,
    read: R,
) -> Result<GenomeBuild, GenomeAssemblyParseError>
where
    R: BufRead,
{
    let mut contigs = vec![];

    for (i, line) in read.lines().enumerate() {
        let line = line?;

        if !line.starts_with("#") {
            let fields: Vec<_> = line.split("\t").collect();
            let record = GenomeAssemblyRecord::try_from(fields)
                .map_err(|e| GenomeAssemblyParseError::AssemblyReportLineParseError(i, e.into()))?;
            contigs.push(Contig::try_from(record).map_err(|e| {
                GenomeAssemblyParseError::AssemblyReportLineParseError(i, e.into())
            })?);
        }
    }

    Ok(GenomeBuild::new(id, contigs))
}

/// Represents the reasons for which parsing of the genome assembly report can fail.
///
/// The reasons include I/O error of the underlying [`BufRead`]
/// (reported as [`GenomeAssemblyParseError::IoError`]) and parsing of the assembly report lines
/// (reported as [`GenomeAssemblyParseError::AssemblyReportLineParseError`]) with one of the following:
/// * Less than 10 columns in the assembly report data line (`Sequence-Name`)
/// * Missing/unparsable column `9` (`Sequence-Length`)
/// * Sequence length being negative (should not really happen)
#[derive(Debug, Error)]
pub enum GenomeAssemblyParseError {
    #[error("I/O error: {0}")]
    IoError(#[from] io::Error),
    #[error("Line #{0}: {1}")]
    AssemblyReportLineParseError(usize, Box<dyn std::error::Error>),
}

struct GenomeAssemblyRecord<'a> {
    fields: Vec<&'a str>,
}

const SEQ_NAME_COL_IDX: usize = 0;
const GENBANK_ACC_COL_IDX: usize = 4;
const REFSEQ_ACC_COL_IDX: usize = 6;
const UCSC_ACC_COL_IDX: usize = 9;
const CONTIG_LEN_COL_IDX: usize = 8;

impl<'a> GenomeAssemblyRecord<'a> {
    fn sequence_name(&self) -> &str {
        self.fields.get(SEQ_NAME_COL_IDX).unwrap()
    }

    fn genbank_acc(&self) -> Option<&str> {
        self.fields
            .get(GENBANK_ACC_COL_IDX)
            .filter(|x| **x != "na")
            .copied()
    }

    fn refseq_acc(&self) -> Option<&str> {
        self.fields
            .get(REFSEQ_ACC_COL_IDX)
            .filter(|x| **x != "na")
            .copied()
    }

    fn ucsc_acc(&self) -> Option<&str> {
        self.fields
            .get(UCSC_ACC_COL_IDX)
            .filter(|x| **x != "na")
            .copied()
    }

    fn length(&self) -> &str {
        self.fields
            .get(CONTIG_LEN_COL_IDX)
            .expect("We checked that there are at least 10 fields in the record!")
    }
}

impl<'a> TryFrom<Vec<&'a str>> for GenomeAssemblyRecord<'a> {
    type Error = &'static str;

    fn try_from(fields: Vec<&'a str>) -> Result<Self, Self::Error> {
        if fields.len() < 10 {
            Err("less than 10 fields in the line")
        } else {
            Ok(GenomeAssemblyRecord { fields })
        }
    }
}

impl<'a> TryFrom<GenomeAssemblyRecord<'a>> for Contig {
    type Error = ParseIntError;
    fn try_from(value: GenomeAssemblyRecord<'a>) -> Result<Self, Self::Error> {
        Ok(Contig {
            name: value.sequence_name().to_string(),
            genbank_name: value.genbank_acc().map(ToString::to_string),
            refseq_name: value.refseq_acc().map(ToString::to_string),
            ucsc_name: value.ucsc_acc().map(ToString::to_string),
            length: value.length().parse()?,
        })
    }
}
