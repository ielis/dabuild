# dabuild

`dabuild` provides your analysis with genome build metadata.

## Example

Use GRCh38.p13 build (*Homo sapiens*):

```rust
use dabuild::{GenomeBuild, GenomeBuildIdentifier};
use dabuild::builds::get_grch38_p13;

// Load the build
let build: GenomeBuild<u32> = get_grch38_p13();

// Check the basic credentials, such as major assembly and patch version
assert_eq!(build.id().major_assembly(), "GRCh38");
assert_eq!(build.id().patch(), Some("p13"));

// Obtain a contig (e.g. `chrY`) by name ...
let y = build.contig_by_name("Y");
assert!(y.is_some());

/// ... or GenBank accession ...
let y = build.contig_by_name("CM000686.2");
assert!(y.is_some());

/// ... or RefSeq accession ...
let y = build.contig_by_name("NC_000024.10");
assert!(y.is_some());

/// ... or UCSC identifier.
let y = build.contig_by_name("chrY");
assert!(y.is_some());
```

## Documentation

See more examples along with the complete documentation at [docs.rs](https://docs.rs/dabuild/latest/dabuild/).
