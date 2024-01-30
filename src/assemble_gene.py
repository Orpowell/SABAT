import tempfile
import pandas as pd
import subprocess
import sys
import os
import warnings
from Bio import SeqIO
from Bio.Seq import Seq

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
        
        self.exon_data["name"] = pd.Categorical(self.exon_data["name"], categories=self.exon_list, ordered=True)
        self.exon_data.sort_values("name", inplace=True)

        self.blastdb: str = BLASTDB_PATH

        self.batch_entry_file = tempfile.NamedTemporaryFile(delete=False)

        self.exon_sequence_file = tempfile.NamedTemporaryFile(delete=False)

        self.protein = ""

        self.cds = ""

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

        if len(self.exon_data.strand.unique()) != 1:
            self.nuke()
            print("Error: all exon strands must be identical!")
            sys.exit(1)

        with open(self.exon_sequence_file.name) as handle:
            self.exon_data["sequence"] = [record.seq for record in SeqIO.parse(handle, "fasta")]

            if self.exon_data.strand.unique()[0] == "-":
                self.exon_data.sequence = self.exon_data.sequence.map(lambda x: x.reverse_complement())

        self.exon_data["ORF1"] = self.exon_data.sequence.map(lambda x: x if len(x) % 3 == 0 else x[:-(len(x) % 3)])
        self.exon_data["ORF2"] = self.exon_data.sequence.map(lambda x: x[1:] if len(x[1:]) % 3 == 0 else x[1:-(len(x[1:]) % 3)])
        self.exon_data["ORF3"] = self.exon_data.sequence.map(lambda x: x[2:] if len(x[2:]) % 3 == 0 else x[2:-(len(x[2:]) % 3)])
    
    def translate_ORFS(self):
        self.exon_data["prot1"] = self.exon_data.ORF1.map(lambda x: x.translate(to_stop=True))
        self.exon_data["prot2"] = self.exon_data.ORF2.map(lambda x: x.translate(to_stop=True))
        self.exon_data["prot3"] = self.exon_data.ORF3.map(lambda x: x.translate(to_stop=True))


    def predict_protein(self):
        protein = []
        cds = []

        for index, exon in self.exon_data.iterrows():
            exon_cds = [exon.ORF1, exon.ORF2, exon.ORF3]
            exon_prots = [exon.prot1, exon.prot2, exon.prot3]
            prot_len = list(map(len, exon_prots))

            biggest_exon = prot_len.index(max(prot_len))
            
            protein.append(exon_prots[biggest_exon])
            cds.append(exon_cds[biggest_exon])
        
        self.protein = "".join([str(seq) for seq in protein])
        self.cds = "".join([str(seq) for seq in cds])

        print(self.protein, "\n")
        print(self.cds)

    def generate_statistics(self) -> None:
        print(f"\nPredicted coverage: {self.exon_data.score.sum()}")
        print(f"Protein length: {len(self.protein)}")
        print(f"CDS gene length: {len(self.cds)}")

    def nuke(self) -> None:
        try:
            os.remove(self.batch_entry_file.name)
        except FileNotFoundError:
            print("batch file not found")

        try:
            os.remove(self.exon_sequence_file.name)
        except FileNotFoundError:
            print("sequences file not found")

        print("Temporary files deleted...")

    def run(self) -> None:

        self.extract_exon_sequences()
        self.load_ORFS_into_dataframe()
        self.translate_ORFS()
        self.predict_protein()
        self.generate_statistics()
        self.nuke()


if __name__ == "__main__":
    CS_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/wheat_genomes/chinese_spring/cs_blastdb/cs_blastdb"
    LONG_BLASTDB_PATH = "/home/powellor/Documents/projects/sr62_homeologues/data/sitopsis_genomes/longissima_blastdb/longissima_blastdb"
    NLR_BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_nlr_cds_refined.mega.bed"
    KINASE_BED_FILE = "/home/powellor/Documents/projects/sr62_homeologues/scripts/SABAT_pipeline/SABAT/src/longissima_sr62_cds_refined.mega.bed"
    NLR_EXON_LIST = ["exon_42", "exon_43", "exon_44"]
    KINASE_EXON_LIST = ["exon_11", "exon_12", "exon_13", "exon_14", "exon_15", "exon_16", "exon_17", "exon_18", "exon_19", "exon_20", "exon_21"]
    KINASE_EXON_LIST = KINASE_EXON_LIST[::-1]
    print(KINASE_EXON_LIST)

    print("\nAnalysing NLR\n")
    NLR_gene = GeneAssembler(
        BED_FILE=NLR_BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=NLR_EXON_LIST
    )
    NLR_gene.run()


    print("\nAnalysing Kinase\n")
    KINASE_gene = GeneAssembler(
        BED_FILE=KINASE_BED_FILE, BLASTDB_PATH=LONG_BLASTDB_PATH, EXON_LIST=KINASE_EXON_LIST
    )
    KINASE_gene.run()
