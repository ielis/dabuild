use std::{error::Error, fs::File, io::BufReader, str::FromStr};

use dabuild::{builds::*, GenomeBuild, GenomeBuildIdentifier};

#[test]
fn grch38_p13() {
    let build = get_grch38_p13::<usize>();

    assert_eq!(build.id().major_assembly(), "GRCh38");
    assert_eq!(build.id().patch(), Some("p13"));
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

#[test]
fn test_parse_assembly_report() -> Result<(), Box<dyn Error>> {
    let path = "data/GCF_000001635.27_GRCm39_assembly_report.txt";
    let read = BufReader::new(File::open(path)?);
    let build = parse_assembly_report(GenomeBuildIdentifier::from_str("GRCm39").unwrap(), read);

    assert!(build.is_ok());
    let build: GenomeBuild<u32> = build?;

    assert_eq!(build.id().major_assembly(), "GRCm39");
    assert_eq!(build.id().patch(), None);
    assert_eq!(build.contigs().len(), 61);

    let contig = build.contig_by_name("Y");
    assert!(contig.is_some());

    let contig = contig.unwrap();
    assert_eq!(contig.name(), "Y");

    let alt: Vec<_> = contig.alt_names().collect();
    assert_eq!(alt.len(), 2);
    assert!(["CM001014.3", "NC_000087.8"]
        .iter()
        .all(|x| alt.contains(x)));

    assert_eq!(contig.length(), &91_455_967u32);

    Ok(())
}
