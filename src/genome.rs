//! # Contigs and genome builds
//!
//! The module includes reference genome build elements such as [`Contig`] and [`GenomeBuild`].

/* ***************************************************************************************************************** *
 *                                               Contig
 * ***************************************************************************************************************** */

use std::str::FromStr;

use num_traits::Zero;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Contig<C> {
    name: String,
    alt_names: Vec<String>,
    length: C,
}

impl<C> Contig<C> {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn alt_names(&self) -> impl Iterator<Item = &str> {
        self.alt_names.iter().map(AsRef::as_ref)
    }

    pub fn length(&self) -> &C {
        &self.length
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

/* ***************************************************************************************************************** *
 *                                               Genome Build
 * ***************************************************************************************************************** */

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomeBuildIdentifier {
    major_assembly: String,
    patch: Option<String>,
}

/// Create identifier from a string.
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
    pub fn major_assembly(&self) -> &str {
        &self.major_assembly
    }
    pub fn patch(&self) -> Option<&str> {
        self.patch.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeBuild<C> {
    id: GenomeBuildIdentifier,
    contigs: Vec<Contig<C>>,
}

impl<C> GenomeBuild<C> {
    pub fn new(id: GenomeBuildIdentifier, mut contigs: Vec<Contig<C>>) -> Self {
        contigs.sort_by(|l, r| l.name().cmp(r.name()));
        GenomeBuild { id, contigs }
    }

    pub fn id(&self) -> &GenomeBuildIdentifier {
        &self.id
    }

    pub fn contigs(&self) -> &[Contig<C>] {
        &self.contigs
    }

    pub fn contig_by_name(&self, name: &str) -> Option<&Contig<C>> {
        self.contigs
            .iter()
            .find(|&c| c.name().eq(name) || c.alt_names().any(|alt_name| alt_name.eq(name)))
    }
}
