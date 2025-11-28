use dabuild::Contig;

#[test]
fn contig_basics() {
    let contig = Contig::builder()
        .name("1")
        .length(248_956_422u32)
        .genbank_name("CM000663.2")
        .refseq_name("NC_000001.11")
        .ucsc_name("chr1")
        .build();

    assert_eq!(contig.name(), "1");
    assert_eq!(contig.alt_names().count(), 3);
    assert_eq!(
        contig.alt_names().collect::<Vec<_>>(),
        vec!["CM000663.2", "NC_000001.11", "chr1"]
    );
    assert_eq!(contig.genbank_name(), Some("CM000663.2"));
    assert_eq!(contig.refseq_name(), Some("NC_000001.11"));
    assert_eq!(contig.ucsc_name(), Some("chr1"));
    assert_eq!(contig.length(), 248_956_422);
}
