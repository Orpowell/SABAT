import tempfile
import pandas as pd
import subprocess
import sys
import os
import warnings
from Bio import SeqIO
import click

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
                "Name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
            ],
        )

        self.exon_list: list[str] = EXON_LIST

        self.exon_data = self.bed[self.bed.Name.isin(self.exon_list)]

        self.exon_data["Name"] = pd.Categorical(
            self.exon_data["Name"].copy(), categories=self.exon_list, ordered=True
        )
        self.exon_data.sort_values("Name", inplace=True)

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
            self.exon_data["sequence"] = [
                record.seq for record in SeqIO.parse(handle, "fasta")
            ]

            if self.exon_data.strand.unique()[0] == "-":
                self.exon_data.sequence = self.exon_data.sequence.map(
                    lambda x: x.reverse_complement()
                )

    def trim_ORFs(self):
        self.exon_data["ORF1"] = self.exon_data.sequence.map(
            lambda x: x if len(x) % 3 == 0 else x[: -(len(x) % 3)]
        )
        self.exon_data["ORF2"] = self.exon_data.sequence.copy().map(
            lambda x: x[1:] if len(x[1:]) % 3 == 0 else x[1 : -(len(x[1:]) % 3)]
        )
        self.exon_data["ORF3"] = self.exon_data.sequence.map(
            lambda x: x[2:] if len(x[2:]) % 3 == 0 else x[2 : -(len(x[2:]) % 3)]
        )

        self.exon_data.ORF1 = self.exon_data.ORF1.map(
            lambda x: x[3:] if str(x).startswith(("TAG", "TAA", "TGA")) else x
        )
        self.exon_data.ORF2 = self.exon_data.ORF2.map(
            lambda x: x[3:] if str(x).startswith(("TAG", "TAA", "TGA")) else x
        )
        self.exon_data.ORF3 = self.exon_data.ORF3.map(
            lambda x: x[3:] if str(x).startswith(("TAG", "TAA", "TGA")) else x
        )

    def translate_ORFS(self):
        self.exon_data["prot1"] = self.exon_data.ORF1.map(
            lambda x: x.translate(to_stop=True)
        )
        self.exon_data["prot2"] = self.exon_data.ORF2.map(
            lambda x: x.translate(to_stop=True)
        )
        self.exon_data["prot3"] = self.exon_data.ORF3.map(
            lambda x: x.translate(to_stop=True)
        )

    def predict_protein(self):
        protein = []
        cds = []
        first_exon = self.exon_list[0]
        last_exon = self.exon_list[-1]

        for index, exon in self.exon_data.iterrows():
            exon_cds = [exon.ORF1, exon.ORF2, exon.ORF3]
            exon_prots = [exon.prot1, exon.prot2, exon.prot3]
            prot_len = list(map(len, exon_prots))

            if exon.Name == first_exon:
                valid_starts = [
                    exon_cds.index(e) for e in exon_cds if str(e).startswith("ATG")
                ]

                if len(valid_starts) == 0:
                    self.nuke()
                    print("No valid start codon in first exon")
                    sys.exit(1)

                elif len(valid_starts) == 1:
                    biggest_exon = valid_starts[0]

            elif exon.Name == last_exon:
                valid_stops = [
                    exon_cds.index(e)
                    for e in exon_cds
                    if str(e).endswith(("TGA", "TAA", "TAG"))
                ]

                if len(valid_stops) == 0:
                    self.nuke()
                    print("No valid stop codon in first exon")
                    sys.exit(1)

                elif len(valid_stops) == 0:
                    biggest_exon = valid_stops[0]

                else:
                    valid_prot_len = [
                        length if n in valid_stops else 0
                        for n, length in enumerate(prot_len)
                    ]
                    biggest_exon = valid_prot_len.index(max(valid_prot_len))

            else:
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
        self.trim_ORFs()
        self.translate_ORFS()
        self.predict_protein()
        self.generate_statistics()
        self.nuke()


@click.command()
@click.option("-i", "--input", type=click.Path(exists=True), required=True, help="bed file")
@click.option("-db", "--blastdb", required=True, help="Path to the BLASTDB (including name)")
@click.option("-e", "--exons", type=click.Path(exists=True), required=True, help="txt file with exons of interest")
def assemble_gene(input: str, blastdb: str, exons: list[str]):
    with open(exons) as file:
        exons_list = [line.strip() for line in file]
        print(exons_list)

    gene = GeneAssembler(BED_FILE=input, BLASTDB_PATH=blastdb, EXON_LIST=exons_list)
    gene.run()
