//! A module with the most frequently used genome builds.
//!
//! We provide these builds out-of-shelf:
//! * GRCh37.p13
//! * GRCh38.p13
//!
//! Other genome builds can be read with [`crate::genome::builds::parse_assembly_report`] function.

use std::{error::Error, io::BufRead, str::FromStr};

use num_traits::Zero;

use super::{Contig, GenomeBuild, GenomeBuildIdentifier};

#[allow(non_upper_case_globals)]
const GRCh37_p13: &[u8] = include_bytes!("data/GCF_000001405.25_GRCh37.p13_assembly_report.tsv");
#[allow(non_upper_case_globals)]
const GRCh38_p13: &[u8] = include_bytes!("data/GCF_000001405.39_GRCh38.p13_assembly_report.tsv");

/// Get the *GRCh37.p13* build.
pub fn get_grch37_p13<C>() -> GenomeBuild<C>
where
    C: FromStr + Zero + PartialOrd,
{
    let id = GenomeBuildIdentifier::new("GRCh37".into(), "p13".into());
    parse_assembly_report(id, GRCh37_p13).expect("Reading builtin GRCh37.p13 assembly report")
}

/// Get the *GRCh38.p13* build.
pub fn get_grch38_p13<C>() -> GenomeBuild<C>
where
    C: FromStr + Zero + PartialOrd,
{
    let id = GenomeBuildIdentifier::new("GRCh38".into(), "p13".into());
    parse_assembly_report(id, GRCh38_p13).expect("Reading builtin GRCh38.p13 assembly report")
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
/// The parsing can fail from several reasons:
///
/// * I/O error of the underlying [`BufRead`]
/// * Missing column `0` (`Sequence-Name`)
/// * Missing/unparsable column `8` (`Sequence-Length`)
/// * Sequence length being negative (should not really happen)
pub fn parse_assembly_report<C, R>(
    id: GenomeBuildIdentifier,
    read: R,
) -> Result<GenomeBuild<C>, Box<dyn Error>>
where
    C: FromStr + Zero + PartialOrd,
    R: BufRead,
{
    let mut contigs = vec![];

    for (i, line) in read.lines().enumerate() {
        // Bail in case of I/O errors.
        let line = line?;

        if line.starts_with("#") {
            continue;
        }
        let fields: Vec<_> = line.split("\t").collect();

        // Disabling the lint to emphasize accessing the columns with indices.
        #[allow(clippy::get_first)]
        let name = if let Some(&name) = fields.get(0) {
            name
        } else {
            return Err(format!("Missing column #0 (`Sequence-Name`) in line #{i} {line}").into());
        };
        let mut alt_names = vec![];

        // Accessions:
        // GenBank, column #4
        if let Some(&gen_bank) = fields.get(4) {
            alt_names.push(gen_bank);
        };
        // RefSeq, column #6
        if let Some(&refseq) = fields.get(6) {
            alt_names.push(refseq);
        };
        // UCSC, column #9
        if let Some(&ucsc) = fields.get(9) {
            alt_names.push(ucsc);
        };

        // Length
        let length = if let Some(&l) = fields.get(8) {
            match l.parse() {
                Ok(length) => length,
                Err(_) => {
                    return Err(format!("Cannot parse field #8 {l:?} into contig length").into())
                }
            }
        } else {
            return Err(
                format!("Missing column #8 (`Sequence-Length`) in line #{i} {line}").into(),
            );
        };

        match Contig::new(name, &alt_names, length) {
            Some(contig) => contigs.push(contig),
            None => return Err("Cannot parse contig".into()),
        };
    }

    Ok(GenomeBuild::new(id, contigs))
}
