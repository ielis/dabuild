//! # Contigs and genome builds
//!
//! The module includes reference genome build elements such as [`Contig`] and [`GenomeBuild`].

/* ***************************************************************************************************************** *
 *                                               Contig
 * ***************************************************************************************************************** */

use std::str::FromStr;

use num_traits::{CheckedSub, Zero};

/// The contig data, such as identifiers and its length.
///
/// `C` is the data type to represent the number of contig's base pairs.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Contig<C> {
    name: String,
    alt_names: Vec<String>,
    length: C,
}

impl<C> Contig<C> {
    /// Get the main name of the contig (e.g. `10`, `X`, `MT`).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the alternative contig identifiers.
    ///
    /// For instance, `CM000686.2`, `NC_000024.10`, and `chrY` for chromosome `Y`.
    pub fn alt_names(&self) -> impl Iterator<Item = &str> {
        self.alt_names.iter().map(AsRef::as_ref)
    }

    /// Get the number of bases of the contig
    pub fn length(&self) -> &C {
        &self.length
    }

    /// Transpose coordinate on a double-stranded sequence to the opposite strand.
    ///
    /// Returns `None` if the operation would lead to underflow.
    pub fn transpose_coordinate(&self, other: &C) -> Option<C>
    where
        C: CheckedSub,
    {
        self.length
            .checked_sub(other)
    }
}

impl<C> Contig<C>
where
    C: Zero + PartialOrd,
{
    pub fn new<T, U>(name: T, alt_names: &[U], length: C) -> Option<Self>
    where
        T: ToString,
        U: ToString,
    {
        if length < C::zero() {
            None
        } else {
            Some(Self {
                name: name.to_string(),
                alt_names: alt_names.iter().map(ToString::to_string).collect(),
                length,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Contig;

    #[test]
    fn test_transpose_coordinate() {
        let contig = Contig::new("X", &["Y"], 10u8).unwrap();

        assert_eq!(contig.transpose_coordinate(&10).unwrap(), 0u8);
        assert_eq!(contig.transpose_coordinate(&8).unwrap(), 2u8);
    }

    #[test]
    fn test_transpose_coordinate_panics() {
        let contig = Contig::new("X", &["Y"], 10u8).unwrap();

        assert!(contig.transpose_coordinate(&11).is_none())
    }
}

/* ***************************************************************************************************************** *
 *                                               Genome Build
 * ***************************************************************************************************************** */

/// Includes information to identify a genome build.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomeBuildIdentifier {
    major_assembly: String,
    patch: Option<String>,
}

/// Create [`GenomeBuildIdentifier`] from a `&str`,
/// using it as a major assembly.
///
/// Infallible.
impl FromStr for GenomeBuildIdentifier {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(GenomeBuildIdentifier {
            major_assembly: s.to_string(),
            patch: None,
        })
    }
}

/// Create [`GenomeBuildIdentifier`] from a tuple.
///
/// The tuple must contain two items:
/// * major assembly
/// * patch
///
/// Use [`GenomeBuildIdentifier::from_str`]
/// to create the identifier without a patch.
impl<T> From<(T, T)> for GenomeBuildIdentifier
where
    T: ToString,
{
    fn from(value: (T, T)) -> Self {
        GenomeBuildIdentifier {
            major_assembly: value.0.to_string(),
            patch: Some(value.1.to_string()),
        }
    }
}

impl GenomeBuildIdentifier {
    /// Get a `&str` with the major assembly identifier.
    pub fn major_assembly(&self) -> &str {
        &self.major_assembly
    }

    /// Get the patch identifier
    /// or `None` if the build identifier has no patch info.
    pub fn patch(&self) -> Option<&str> {
        self.patch.as_deref()
    }
}

/// Genome build includes the contigs and genome build metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeBuild<C> {
    id: GenomeBuildIdentifier,
    contigs: Vec<Contig<C>>,
}

impl<C> GenomeBuild<C> {
    pub fn new<I>(id: GenomeBuildIdentifier, contigs: I) -> Self
    where
        I: IntoIterator<Item = Contig<C>>,
    {
        let mut contigs: Vec<_> = contigs.into_iter().collect();
        contigs.sort_by(|l, r| l.name().cmp(r.name()));
        GenomeBuild { id, contigs }
    }

    /// Get the genome build identifiers.
    pub fn id(&self) -> &GenomeBuildIdentifier {
        &self.id
    }

    /// Get an iterator with all contigs.
    pub fn contigs(&self) -> impl Iterator<Item = &Contig<C>> {
        self.contigs.iter()
    }

    pub fn contig_by_name(&self, name: &str) -> Option<&Contig<C>> {
        self.contigs
            .iter()
            .find(|&c| c.name().eq(name) || c.alt_names().any(|alt_name| alt_name.eq(name)))
    }
}
