import pandas as pd


def swap_sstart_send(df):
    cond = df.sstart > df.send
    df.loc[cond, ["sstart", "send"]] = df.loc[cond, ["send", "sstart"]].values
    return df


def BLAST2BED9(input, output):
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

    # Establish oreintation of BLAST hits
    cds_blast_data["orientation"] = cds_blast_data.sstart < cds_blast_data.send

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
    if any(cds_blast_data.sseqid.str.contains("|")):
        bed9["chrom"] = cds_blast_data.sseqid.map(lambda x: x.split("|")[1])

    else:
        bed9["chrom"]: str = cds_blast_data.sseqid

    # Fill columns 2-9
    bed9["chromStart"]: int = cds_blast_data.sstart
    bed9["chromEnd"]: int = cds_blast_data.send 
    bed9["name"]: int = cds_blast_data.index
    bed9["score"]: float = round(cds_blast_data.length / cds_blast_data.qlen, 2)
    bed9["strand"]: str = cds_blast_data.orientation.map(
        lambda x: "+" if x is True else "-"
    )
    bed9["thickStart"]: int = cds_blast_data.sstart 
    bed9["thickEnd"]: int = cds_blast_data.send
    bed9["itemRgb"]: str = "145,30,180"

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
    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds.blastn",
        "Sr62_cds.bed",
    )
    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds_e0001.blastn",
        "Sr62_cds_refined.bed",
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/SrNLR_cds.blastn",
        "SrNLR_cds.bed",
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/NLR_cds_e0001.blastn",
        "SrNLR_cds_refined.bed",
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/NLR_cds_e0001.blastn.dcmega",
        "SrNLR_cds_refined.mega.bed",
    )

    BLAST2BED9(
        "~/Documents/projects/sr62_homeologues/analysis/blast_annotation_mapping/test_data/chinese_spring/Sr62_cds_e0001.blastn.dcmega",
        "Sr62_cds_refined.mega.bed",
    )
