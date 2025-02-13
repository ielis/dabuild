use dabuild::Contig;

#[test]
fn contig_basics() {
    let contig = Contig::new("1", &["CM000663.2", "NC_000001.11", "chr1"], 10u8).unwrap();

    assert_eq!(contig.name(), "1");
    assert_eq!(contig.alt_names().count(), 3);
    assert_eq!(
        contig.alt_names().collect::<Vec<_>>(),
        vec!["CM000663.2", "NC_000001.11", "chr1"]
    );
    assert_eq!(contig.length(), &10u8);
}
