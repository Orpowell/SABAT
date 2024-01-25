import os
import pandas as pd
from Bio import SeqIO


def fetch_exons(BED_FILE, EXON_LIST, BLASTDB_PATH) -> None:
    bed = pd.read_csv(
        BED_FILE,
        sep="\t",
        header=None,
        names=[
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

    bed = bed[bed.name.isin(EXON_LIST)]

    with open("tmp_entry_batch.txt", "w+") as file:
        for index, row in bed.iterrows():
            file.write(f"{row.chrom} {row.chromStart}-{row.chromEnd}\n")

    blastcmd = f" blastdbcmd -db {BLASTDB_PATH} -entry_batch tmp_entry_batch.txt -out exon_sequences.fa -outfmt %f"

    os.system(blastcmd)


def predict_gene() -> None:
    with open("exon_sequences.fa") as handle:
        seqs = [record.seq for record in SeqIO.parse(handle, "fasta")]

    predicted_exons = []

    for n, seq in enumerate(seqs):
        sequence = ""
        exon_length = 0

        for i in range(0, 3):
            exon = seq[i:].translate(to_stop=True)

            if len(exon) > exon_length:
                sequence = exon
                exon_length = len(exon)

        predicted_exons.append(sequence)

    print("".join([str(seq) for seq in predicted_exons]))


if __name__ == "__main__":
    BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/wheat_genomes/chinese_spring/cs_blastdb/cs_blastdb"

    # simple nlr BLAST
    BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/SrNLR_cds.bed"
    EXON_LIST = [6, 7, 8]
    fetch_exons(BED_FILE=BED_FILE, EXON_LIST=EXON_LIST, BLASTDB_PATH=BLASTDB_PATH)
    predict_gene()

    # refined NLR
    BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/SrNLR_cds_refined.bed"
    EXON_LIST = [118, 119, 120]
    fetch_exons(BED_FILE=BED_FILE, EXON_LIST=EXON_LIST, BLASTDB_PATH=BLASTDB_PATH)
    predict_gene()
