//! # Contigs and genome builds
//!
//! The module includes reference genome build elements such as [`Contig`] and [`GenomeBuild`].

/* ***************************************************************************************************************** *
 *                                               Contig
 * ***************************************************************************************************************** */
use std::convert::Infallible;
use std::str::FromStr;

/// The contig data, such as identifiers and its length.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Contig {
    pub(crate) name: String,
    pub(crate) genbank_name: Option<String>,
    pub(crate) refseq_name: Option<String>,
    pub(crate) ucsc_name: Option<String>,
    pub(crate) length: u32,
}

impl Contig {
    /// Get the main name of the contig (e.g. `10`, `X`, `MT`).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the alternative contig identifiers.
    ///
    /// For instance, `CM000686.2`, `NC_000024.10`, and `chrY` for chromosome `Y`.
    ///
    /// # Example
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["CM000686.2", "NC_000024.10", "chrY"], 57_227_415).expect("The contig data are valid");
    ///
    /// let alt_names: Vec<_> = contig.alt_names().collect();
    /// assert_eq!(&alt_names, &["CM000686.2", "NC_000024.10", "chrY"]);
    pub fn alt_names(&self) -> impl Iterator<Item = &str> {
        std::iter::once(&self.genbank_name)
            .chain(std::iter::once(&self.refseq_name))
            .chain(std::iter::once(&self.ucsc_name))
            .flatten()
            .map(AsRef::as_ref)
    }

    /// Get the GenBank contig identifier, if available.
    ///
    /// For instance, `CM000686.2` for chromosome `Y` of the *GRCh38.p13* assembly.
    ///
    /// # Example
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["CM000686.2", "NC_000024.10", "chrY"], 57_227_415).expect("The contig data are valid");
    ///
    /// assert_eq!(contig.genbank_name(), Some("CM000686.2"));
    /// ```
    pub fn genbank_name(&self) -> Option<&str> {
        self.genbank_name.as_deref()
    }

    /// Get the RefSeq contig identifier, if available.
    ///
    /// For instance, `NC_000024.10` for chromosome `Y` of the *GRCh38.p13* assembly.
    ///
    /// # Example
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["CM000686.2", "NC_000024.10", "chrY"], 57_227_415).expect("The contig data are valid");
    ///
    /// assert_eq!(contig.refseq_name(), Some("NC_000024.10"));
    /// ```
    pub fn refseq_name(&self) -> Option<&str> {
        self.refseq_name.as_deref()
    }

    /// Get the UCSC contig identifier, if available.
    ///
    /// For instance, `chrY` for chromosome `Y` of the *GRCh38.p13* assembly.
    ///
    /// # Example
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["CM000686.2", "NC_000024.10", "chrY"], 57_227_415).expect("The contig data are valid");
    ///
    /// assert_eq!(contig.ucsc_name(), Some("chrY"));
    /// ```
    pub fn ucsc_name(&self) -> Option<&str> {
        self.ucsc_name.as_deref()
    }

    /// Get the number of bases of the contig, a.k.a. its length.
    pub fn length(&self) -> u32 {
        self.length
    }

    /// Transpose coordinate on a double-stranded sequence to the opposite strand.
    ///
    /// Returns `None` if the operation would lead to an overflow.
    pub fn transpose_coordinate(&self, other: u32) -> Option<u32> {
        self.length.checked_sub(other)
    }

    /// Create a Contig from `name`, the alternative names, and its length.
    ///
    /// # Example
    ///
    /// Create a `Contig` with fields from the GRC assembly report:
    ///
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["CM000686.2", "NC_000024.10", "chrY"], 57_227_415);
    ///
    /// assert!(contig.is_some());
    /// let contig = contig.unwrap();
    /// assert_eq!(contig.name(), "Y");
    /// assert_eq!(contig.genbank_name(), Some("CM000686.2"));
    /// ```
    ///
    /// Create a `Contig` with sequence name and UCSC accession:
    ///
    /// ```
    /// use dabuild::Contig;
    ///
    /// let contig = Contig::new("Y", &["na", "", "chrY"], 57_227_415);
    ///
    /// assert!(contig.is_some());
    ///
    /// let contig = contig.unwrap();
    /// assert_eq!(contig.name(), "Y");
    /// assert!(contig.genbank_name().is_none());
    /// assert!(contig.refseq_name().is_none());
    /// assert_eq!(contig.ucsc_name(), Some("chrY"));
    /// ```
    ///
    /// The `new` expects `alt_names` with three items that correspond to:
    /// * GenBank accession
    /// * RefSeq accession
    /// * UCSC accession
    ///
    /// An accession equaling to an empty string or `"na"` is filtered out.
    pub fn new(name: impl ToString, alt_names: &[impl ToString], length: u32) -> Option<Self> {
        const NON_EMPTY_NON_NA_STRING: fn(&String) -> bool = |v| !v.is_empty() && v != "na";

        Some(Self {
            name: name.to_string(),
            #[allow(clippy::get_first)]
            genbank_name: alt_names
                .get(0)
                .map(ToString::to_string)
                .filter(NON_EMPTY_NON_NA_STRING),
            refseq_name: alt_names
                .get(1)
                .map(ToString::to_string)
                .filter(NON_EMPTY_NON_NA_STRING),
            ucsc_name: alt_names
                .get(2)
                .map(ToString::to_string)
                .filter(NON_EMPTY_NON_NA_STRING),
            length,
        })
    }
}

#[cfg(test)]
mod contig_tests {
    use super::Contig;

    #[test]
    fn test_transpose_coordinate() {
        let contig = Contig::new("X", &["Y"], 10).unwrap();

        assert_eq!(contig.transpose_coordinate(10).unwrap(), 0);
        assert_eq!(contig.transpose_coordinate(8).unwrap(), 2);
    }

    #[test]
    fn test_transpose_coordinate_panics() {
        let contig = Contig::new("X", &["Y"], 10).unwrap();

        assert!(contig.transpose_coordinate(11).is_none())
    }
}

/* ***************************************************************************************************************** *
 *                                               Genome Build
 * ***************************************************************************************************************** */

/// All information needed to identify a genome build.
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
    type Err = Infallible;

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

/// Genome build includes the contigs and the genome build metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeBuild {
    id: GenomeBuildIdentifier,
    contigs: Vec<Contig>,
}

impl GenomeBuild {
    pub fn new<I>(id: GenomeBuildIdentifier, contigs: I) -> Self
    where
        I: IntoIterator<Item = Contig>,
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
    pub fn contigs(&self) -> impl Iterator<Item = &Contig> {
        self.contigs.iter()
    }

    /// Retrieve a [`Contig`] by its name or [`None`] if no such [`Contig`] can be found.
    ///
    /// # Examples
    ///
    /// ```
    /// use dabuild::{GenomeBuild, GenomeBuildIdentifier};
    /// use dabuild::builds::get_grch38_p13;
    ///
    /// let build: GenomeBuild = get_grch38_p13();
    ///
    /// // Query by a contig name ...
    /// let chrY = build.contig_by_name("Y").unwrap();
    /// assert_eq!(chrY.name(), "Y");
    ///
    /// // ... or by GenBank accession ...
    /// let chrY = build.contig_by_name("CM000686.2").unwrap();
    /// assert_eq!(chrY.genbank_name(), Some("CM000686.2"));
    ///
    /// // ... or by RefSeq accession ...
    /// let chrY = build.contig_by_name("NC_000024.10").unwrap();
    /// assert_eq!(chrY.refseq_name(), Some("NC_000024.10"));
    ///
    /// // ... or by UCSC accession.
    /// let chrY = build.contig_by_name("chrY").unwrap();
    /// assert_eq!(chrY.ucsc_name(), Some("chrY"));
    /// ```
    pub fn contig_by_name(&self, name: &str) -> Option<&Contig> {
        self.contigs
            .iter()
            .find(|&c| c.name().eq(name) || c.alt_names().any(|alt_name| alt_name.eq(name)))
    }
}
