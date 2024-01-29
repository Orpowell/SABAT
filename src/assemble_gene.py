import tempfile
import pandas as pd
import subprocess
import sys
import os
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

    def load_ORFS_into_dataframe(self) -> None:

        with open(self.exon_sequence_file.name) as handle:
            self.exon_data["sequence"] = [str(record.seq) for record in SeqIO.parse(handle, "fasta")]

        self.exon_data["ORF1"] = self.exon_data.sequence.map(lambda x: x if len(x) % 3 == 0 else x[:-(len(x) % 3)])
        self.exon_data["ORF2"] = self.exon_data.sequence.map(lambda x: x[1:] if len(x[1:]) % 3 == 0 else x[1:-(len(x[1:]) % 3)])
        self.exon_data["ORF3"] = self.exon_data.sequence.map(lambda x: x[2:] if len(x[2:]) % 3 == 0 else x[2:-(len(x[2:]) % 3)])

        print(self.exon_data)


    def predict_protein(self) -> None:

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

    def nuke(self) -> None:
        os.remove(self.batch_entry_file.name)
        os.remove(self.exon_sequence_file.name)
        print("Temporary files deleted...")

    def run(self) -> None:

        self.extract_exon_sequences()
        #self.predict_protein()
        self.load_ORFS_into_dataframe()
        self.generate_statistics()
        self.nuke()


if __name__ == "__main__":
    CS_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/wheat_genomes/chinese_spring/cs_blastdb/cs_blastdb"
    LONG_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/sitopsis_genomes/longissima_blastdb/longissima_blastdb"
    NLR_BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_nlr_cds_refined.mega.bed"
    KINASE_BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_sr62_cds_refined.mega.bed"
    NLR_EXON_LIST = ["exon_42", "exon_43", "exon_44"]
    KINASE_EXON_LIST = ["exon_11", "exon_12", "exon_13", "exon_14", "exon_15", "exon_16", "exon_17", "exon_18", "exon_19", "exon_20", "exon_21"]

    print("\nAnalysing NLR\n")
    NLR_gene = GeneAssembler(
        BED_FILE=NLR_BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=NLR_EXON_LIST
    )
    NLR_gene.run()

    """
    print("\nAnalysing Kinase\n")
    KINASE_gene = GeneAssembler(
        BED_FILE=KINASE_BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=KINASE_EXON_LIST
    )
    KINASE_gene.run()
    """
