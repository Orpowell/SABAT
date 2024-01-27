import pandas as pd


def BLAST2BED9(
    input: str, output: str, locus_size: int, exon_count: int, q_cov_threshold: float, refseq=False
) -> None:
    cds_blast_data = pd.read_csv(
        input,
        sep="\t",
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
        ],
    )

    # Establish strand orientation and query coverage of BLAST hits
    cds_blast_data["orientation"] = cds_blast_data.sstart < cds_blast_data.send
    cds_blast_data["strand"] = cds_blast_data.orientation.map(
        lambda x: "+" if x is True else "-"
    )
    cds_blast_data["qcov"] = round(cds_blast_data.length / cds_blast_data.qlen, 2)

    # Correctly order start/end of BLAST hits for BED
    cond = cds_blast_data.sstart > cds_blast_data.send
    cds_blast_data.loc[cond, ["sstart", "send"]] = cds_blast_data.loc[
        cond, ["send", "sstart"]
    ].values

    # Sort values + reindex
    cds_blast_data.sort_values(by="sstart", inplace=True)
    cds_blast_data.reset_index(drop=True, inplace=True)

    # Initialise new dataframe for bed file
    bed9 = pd.DataFrame()

    # Remove any formatting from BLASTDB from sequence IDs
    if refseq:
        cds_blast_data["chromosome"]: str = cds_blast_data.sseqid.map(   # type: ignore
            lambda x: x.split("|")[1]
        )

    else:
        cds_blast_data["chromosome"] = cds_blast_data.sseqid

    # Fill columns 2-9
    bed9["chrom"]: str = cds_blast_data.chromosome  # type: ignore
    bed9["chromStart"]: int = cds_blast_data.sstart  # type: ignore
    bed9["chromEnd"]: int = cds_blast_data.send  # type: ignore
    bed9["name"]: str = [f"exon_{i}" for i in cds_blast_data.index]  # type: ignore
    bed9["score"]: float = cds_blast_data.qcov  # type: ignore
    bed9["strand"]: str = cds_blast_data.strand  # type: ignore
    bed9["thickStart"]: int = cds_blast_data.sstart  # type: ignore
    bed9["thickEnd"]: int = cds_blast_data.send  # type: ignore
    bed9["itemRgb"]: str = "145,30,180"  # type: ignore

    # Label gene loci according to parameters and add to bed file
    # This provides confidence in the absence of annotations
    gene_locus = 0

    for window in cds_blast_data.sort_values(["chromosome", "strand"]).rolling(
        exon_count
    ):
        if (
            len(window) > exon_count - 1
            and all(window.chromosome.unique())  # type: ignore
            and all(window.strand.unique())
        ):
            locus_qcov = round(window.qcov.sum(), 2)
            size = window.send.max() - window.sstart.min()

            if locus_qcov > q_cov_threshold and size < locus_size * 1.5:
                bed9.loc[len(bed9.index)] = [
                    window.chromosome.unique()[0],
                    window.sstart.min(),
                    window.send.max(),
                    f"locus_{gene_locus}",
                    locus_qcov,
                    window.strand.unique()[0],
                    window.sstart.min(),
                    window.send.max(),
                    "0,255,0",
                ]
                gene_locus += 1

    # Save output to bed file (formatted as tsv)
    bed9.to_csv(
        output,
        sep="\t",
        header=False,
        index=False,
        columns=[
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
        ],
    )

if __name__ == "__main__":
    NLR_EXONS = 3
    NLR_LOCUS = 6000
    NLR_QCOV = 0.8

    KINASE_EXONS = 9
    KINASE_LOCUS = 10000
    KINASE_QCOV = 0.8

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds.blastn",
        "Sr62_cds.bed",
        exon_count=KINASE_EXONS,
        q_cov_threshold=KINASE_QCOV,
        locus_size=KINASE_LOCUS,
        refseq=True
    )
    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds_e0001.blastn",
        "Sr62_cds_refined.bed",
        exon_count=KINASE_EXONS,
        q_cov_threshold=KINASE_QCOV,
        locus_size=KINASE_LOCUS,
        refseq=True
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/SrNLR_cds.blastn",
        "SrNLR_cds.bed",
        exon_count=NLR_EXONS,
        locus_size=NLR_LOCUS,
        q_cov_threshold=NLR_QCOV,
        refseq=True
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/NLR_cds_e0001.blastn",
        "SrNLR_cds_refined.bed",
        exon_count=NLR_EXONS,
        locus_size=NLR_LOCUS,
        q_cov_threshold=NLR_QCOV,
        refseq=True
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/NLR_cds_e0001.blastn.dcmega",
        "SrNLR_cds_refined.mega.bed",
        exon_count=NLR_EXONS,
        locus_size=NLR_LOCUS,
        q_cov_threshold=NLR_QCOV,
        refseq=True
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds_e0001.blastn.dcmega",
        "Sr62_cds_refined.mega.bed",
        exon_count=KINASE_EXONS,
        q_cov_threshold=KINASE_QCOV,
        locus_size=KINASE_LOCUS,
        refseq=True
    )

    BLAST2BED9(
        "/home/powellor/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/longissima/NLR_cds_e0001.blastn.dcmega",
        "longissima_nlr_cds_refined.mega.bed",
        exon_count=NLR_EXONS,
        q_cov_threshold=NLR_QCOV,
        locus_size=NLR_LOCUS
    )

    BLAST2BED9(
        "/home/powellor/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/longissima/Sr62_cds_e0001.blastn.dcmega",
        "longissima_sr62_cds_refined.mega.bed",
        exon_count=KINASE_EXONS,
        q_cov_threshold=KINASE_QCOV,
        locus_size=KINASE_LOCUS
    )
