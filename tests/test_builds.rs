use dabuild::builds::*;

#[test]
fn grch38_p13() {
    let build = get_grch38_p13::<usize>();

    assert_eq!(build.id().major_assembly(), "GRCh38");
    assert_eq!(build.id().patch(), "p13");
    assert_eq!(build.contigs().len(), 640);

    let contig = build.contig_by_name("chr1");
    assert!(contig.is_some());

    let contig = contig.unwrap();
    assert_eq!(contig.name(), "1");

    let alt: Vec<_> = contig.alt_names().collect();
    assert_eq!(alt.len(), 3);
    assert!(["CM000663.2", "NC_000001.11", "chr1"]
        .iter()
        .all(|x| alt.contains(x)));

    assert_eq!(contig.length(), &248_956_422usize);
}
