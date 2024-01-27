import tempfile
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO


class GeneAssembler:
    def __init__(self, BED_FILE, BLASTDB_PATH, EXON_LIST) -> None:
        self.bed = pd.read_csv(
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

        self.exon_list: list[str] = EXON_LIST

        self.exon_data = self.bed[self.bed.name.isin(self.exon_list)]

        self.blastdb: str = BLASTDB_PATH

        self.sequences_extracted = False

        self.batch_entry_file = tempfile.NamedTemporaryFile(delete=False)

        self.exon_sequence_file = tempfile.NamedTemporaryFile(delete=False)

        print(self.exon_sequence_file.name)

    def extract_exon_sequences(self) -> None:
        print(self.bed)
        print(self.exon_data)
        with open(self.batch_entry_file.name, "w+") as file:
            for index, row in self.exon_data.iterrows():
                file.write(f"{row.chrom} {row.chromStart}-{row.chromEnd}\n")

        try:
            
            subprocess.run(["blastdbcmd", "-db", f"{self.blastdb}", "-entry_batch", f"{self.batch_entry_file.name}", "-out", f"{self.exon_sequence_file.name}", "-outfmt", "%f"], check=True)
            self.sequences_extracted = True

        except (subprocess.CalledProcessError, FileNotFoundError) as error:
            print(error)
            sys.exit(1)

    


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

    blastcmd = f"blastdbcmd -db {BLASTDB_PATH} -entry_batch tmp_entry_batch.txt -out exon_sequences.fa -outfmt %f"

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
    CS_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/wheat_genomes/chinese_spring/cs_blastdb/cs_blastdb"
    LONG_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/sitopsis_genomes/longissima_blastdb/longissima_blastdb"
    BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_nlr_cds_refined.mega.bed"
    EXON_LIST = ["exon_42", "exon_43", "exon_44"]

    fetch_exons(BED_FILE=BED_FILE, EXON_LIST=EXON_LIST, BLASTDB_PATH=LONG_BLASTDB_PATH)
    predict_gene()

    gene = GeneAssembler(BED_FILE=BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=EXON_LIST)
    gene.extract_exon_sequences()