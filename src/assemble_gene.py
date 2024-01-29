import tempfile
import pandas as pd
import subprocess
import sys
import warnings
from Bio import SeqIO

warnings.filterwarnings("ignore")

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

        self.protein = ""

    def extract_exon_sequences(self) -> None:
        with open(self.batch_entry_file.name, "w+") as file:
            for index, row in self.exon_data.iterrows():
                file.write(f"{row.chrom} {row.chromStart}-{row.chromEnd}\n")

        try:
            subprocess.run(
                [
                    "blastdbcmd",
                    "-db",
                    f"{self.blastdb}",
                    "-entry_batch",
                    f"{self.batch_entry_file.name}",
                    "-out",
                    f"{self.exon_sequence_file.name}",
                    "-outfmt",
                    "%f",
                ],
                check=True,
            )
            self.sequences_extracted = True

        except (subprocess.CalledProcessError, FileNotFoundError) as error:
            print(error)
            sys.exit(1)

    def predict_protein(self) -> None:
        if self.sequences_extracted is False:
            print("Error: Sequences not extracted")
            sys.exit(1)

        print(f"Analysing gene containing {len(self.exon_data)} exon(s)\n")

        with open(self.exon_sequence_file.name) as handle:
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

        self.protein = "".join([str(seq) for seq in predicted_exons])

        print(self.protein)

    def generate_statistics(self) -> None:
        print(f"\nPredicted coverage: {self.exon_data.score.sum()}")
        print(f"Protein length: {len(self.protein)}")
        print("CDS gene length: TBA")

    def run(self) -> None:

        self.extract_exon_sequences()
        self.predict_protein()
        self.generate_statistics()


if __name__ == "__main__":
    CS_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/wheat_genomes/chinese_spring/cs_blastdb/cs_blastdb"
    LONG_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/sitopsis_genomes/longissima_blastdb/longissima_blastdb"
    BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_nlr_cds_refined.mega.bed"
    EXON_LIST = ["exon_42", "exon_43", "exon_44"]

    """
    fetch_exons(BED_FILE=BED_FILE, EXON_LIST=EXON_LIST, BLASTDB_PATH=LONG_BLASTDB_PATH)
    predict_gene()
    """

    gene = GeneAssembler(
        BED_FILE=BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=EXON_LIST
    )
    gene.run()
