import tempfile
import pandas as pd
import subprocess
import sys
import os
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from abc import ABC

warnings.filterwarnings("ignore")


class AbstractGeneAssembler(ABC):
    def __init__(self, BED_FILE, BLASTDB_PATH, output, flank) -> None:
        self.input_file = BED_FILE

        self.bed = pd.read_csv(
            self.input_file,
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
                "exon_list",
            ],
        )

        self.exon_data = pd.DataFrame()

        self.blastdb: str = BLASTDB_PATH

        self.output: str = output

        self.batch_entry_file = tempfile.NamedTemporaryFile(delete=False)

        self.exon_sequence_file = tempfile.NamedTemporaryFile(delete=False)

        self.protein = ""

        self.cds = ""

        self.flank = flank

        logging.info(f"Predicting gene from {self.input_file}...")
        logging.info(
            "temporary files created will not be deleted if the program is stopped prematurely..."
        )

        self.strand = ""

    def load_exons_into_dataframe(self) -> None:
        with open(self.exon_sequence_file.name) as handle:
            if self.strand == "-":
                self.exon_data["sequence"] = [
                    record.seq.reverse_complement()
                    for record in SeqIO.parse(handle, "fasta")
                ]
            else:
                self.exon_data["sequence"] = [
                    record.seq for record in SeqIO.parse(handle, "fasta")
                ]

    def extract_exon_sequences(self) -> None:
        logging.info("Generating batch entry file...")

        with open(self.batch_entry_file.name, "w+") as file:
            for index, row in self.exon_data.iterrows():
                file.write(f"{row.chrom} {row.chromStart}-{row.chromEnd}\n")

        logging.info(f"extracting sequences from {self.blastdb}...")

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
            logging.error(error)
            sys.exit(1)

    def filter_exon_data(self) -> None:
        pass

    def extend_flank3(self) -> None:
        if self.strand == "+":
            last_exon = int(self.exon_list[-1].split("_")[1])
            self.exon_data.at[last_exon, "chromEnd"] = (
                self.flank + self.exon_data.at[last_exon, "chromEnd"]
            )

        else:
            last_exon = int(self.exon_list[-1].split("_")[1])
            self.exon_data.at[last_exon, "chromStart"] = (
                self.exon_data.at[last_exon, "chromStart"] - self.flank
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

    def predict_CDS(self):
        first_exon = self.exon_list[0]
        last_exon = self.exon_list[-1]
        cds = []

        for index, exon in self.exon_data.iterrows():
            exon_cds = [exon.ORF1, exon.ORF2, exon.ORF3]
            best_orfs = []
            logging.info(f"analyzing {exon.Name}")
            if exon.Name == first_exon:
                for exon in exon_cds:
                    codons = [exon[i : i + 3] for i in range(0, len(exon), 3)]
                    stops = []

                    try:
                        atg = codons.index("ATG")
                    except ValueError:
                        best_orfs.append("")
                        continue

                    stops.append(len(codons))
                    for terminator in ["TAA", "TGA", "TAG"]:
                        try:
                            stop_pos = codons[atg:].index(terminator)
                            stops.append(stop_pos)
                        except ValueError:
                            continue

                    earliest_stop = min(stops)

                    exon = "".join([str(codon) for codon in codons[atg:earliest_stop]])

                    best_orfs.append(exon)

                lengths = list(map(len, best_orfs))
                max_len = lengths.index(max(lengths))
                exon = best_orfs[max_len]
                cds.append(exon)
                logging.info(exon)

            elif exon.Name == last_exon:
                for exon in exon_cds:
                    codons = [exon[i : i + 3] for i in range(0, len(exon), 3)]
                    ranges = [
                        n
                        for n, codon in enumerate(codons)
                        if codon in ["TAG", "TAA", "TGA"]
                    ]

                    if len(ranges) == 0:
                        best_orfs.append(exon)
                        continue

                    if ranges[-1] != len(codons) - 1:
                        ranges.append(len(codons))

                    start = 0
                    end = 0
                    biggest = []

                    for pos in ranges:
                        end = pos
                        exon = codons[start : end + 1]
                        if len(exon) > len(biggest):
                            biggest = exon

                        start = end + 1

                    best_orfs.append("".join(map(str, biggest)))

                lengths = [
                    len(exon) if str(exon).endswith(("TAG", "TGA", "TAA")) else 0
                    for exon in best_orfs
                ]

                max_len = lengths.index(max(lengths))
                exon = best_orfs[max_len]
                logging.info(exon)
                cds.append(exon)

            else:
                for exon in exon_cds:
                    codons = [exon[i : i + 3] for i in range(0, len(exon), 3)]
                    ranges = [
                        n
                        for n, codon in enumerate(codons)
                        if codon in ["TAG", "TAA", "TGA"]
                    ]

                    if len(ranges) == 0:
                        best_orfs.append(exon)
                        continue

                    if ranges[-1] != len(codons) - 1:
                        ranges.append(len(codons))

                    start = 0
                    end = 0
                    biggest = []

                    for pos in ranges:
                        end = pos
                        exon = codons[start:end]
                        if len(exon) > len(biggest):
                            biggest = exon

                        start = end + 1

                    best_orfs.append("".join(map(str, biggest)))

                lengths = [len(exon) for exon in best_orfs]

                max_len = lengths.index(max(lengths))
                exon = best_orfs[max_len]
                logging.info(exon)
                p = Seq(str(exon)).translate()
                logging.info(p)
                cds.append(exon)

        self.cds = "".join([str(exon) for exon in cds])

    def check_CDS(self):
        if not self.cds.startswith("ATG"):
            self.nuke()

        if not self.cds.endswith(("TAA", "TAG", "TGA")):
            logging.info("CDS doesn't contain a valid stop codon...")
            logging.info("Shutting down")
            self.nuke()
            sys.exit(1)

    def predict_protein(self):
        self.protein = str(Seq(self.cds).translate())

    def write_output_sequences(self):
        logging.info(f"Writing CDS sequence to {self.output}.cds.fasta")
        with open(f"{self.output}.cds.fasta", "w+") as cds:
            cds.write(f">{os.path.basename(self.output)}\n")
            cds.write(self.cds)

        logging.info(f"Writing protein sequence to {self.output}.prot.fasta")
        with open(f"{self.output}.prot.fasta", "w+") as prot:
            prot.write(f">{os.path.basename(self.output)}\n")
            prot.write(self.protein)

    def generate_statistics(self) -> None:
        logging.info(f"Predicted coverage: {self.exon_data.score.sum()}")
        logging.info(f"Protein length: {len(self.protein)}")
        logging.info(f"CDS gene length: {len(self.cds)}")

    def nuke(self) -> None:
        try:
            os.remove(self.batch_entry_file.name)
        except FileNotFoundError:
            logging.error("batch file not found")

        try:
            os.remove(self.exon_sequence_file.name)
        except FileNotFoundError:
            logging.error("sequences file not found")

        logging.info("Temporary files deleted...")

    def run(self) -> None:
        self.filter_exon_data()
        self.extend_flank3()
        self.extract_exon_sequences()
        self.load_exons_into_dataframe()
        self.trim_ORFs()
        self.predict_CDS()
        self.check_CDS()
        self.predict_protein()
        self.write_output_sequences()
        self.generate_statistics()
        self.nuke()


class ExonAssembler(AbstractGeneAssembler):
    def __init__(self, BED_FILE, BLASTDB_PATH, EXON_LIST, output, flank) -> None:
        super().__init__(BED_FILE, BLASTDB_PATH, output, flank)
        self.exon_list: list[str] = EXON_LIST
        logging.info(f"Analysing exons: {' '.join([exon for exon in self.exon_list])}")

    def filter_exon_data(self) -> None:
        exon_data = self.bed[self.bed.Name.isin(self.exon_list)]

        exon_data["Name"] = pd.Categorical(
            exon_data["Name"].copy(), categories=self.exon_list, ordered=True
        )

        exon_data.sort_values("Name", inplace=True)

        if len(exon_data.strand.unique()) != 1:
            self.nuke()
            logging.error("all exon strands must be identical!")
            sys.exit(1)

        self.strand = exon_data.strand.unique()[0]
        self.exon_data = exon_data
        self.strand = self.exon_data.strand.unique()[0]


class LocusAssembler(AbstractGeneAssembler):
    def __init__(self, BED_FILE, BLASTDB_PATH, locus, output, flank) -> None:
        super().__init__(BED_FILE, BLASTDB_PATH, output, flank)
        self.locus = locus

    def filter_exon_data(self) -> None:
        local_min = self.bed.loc[self.bed.Name == self.locus].chromStart.iloc[0]
        local_max = self.bed.loc[self.bed.Name == self.locus].chromEnd.iloc[0]

        exon_data = self.bed[
            (self.bed.chromStart >= local_min)
            & (self.bed.chromEnd <= local_max)
            & (self.bed.Name.map(lambda x: x.startswith("exon")))
        ]

        self.strand = exon_data.strand.unique()[0]

        if self.strand == "-":
            self.exon_list = list(exon_data.Name)[::-1]
            exon_data["Name"] = pd.Categorical(
                exon_data["Name"].copy(), categories=self.exon_list, ordered=True
            )
            exon_data.sort_values(by="Name", inplace=True)

        else:
            self.exon_list = list(exon_data.Name)

        self.exon_data = exon_data
        self.strand = self.exon_data.strand.unique()[0]


def assemble_exons(input: str, blastdb: str, exons: list[str], output: str, flank: int):
    """
    Assemble a gene from exons defined in a bed file
    """
    gene = ExonAssembler(
        BED_FILE=input,
        BLASTDB_PATH=blastdb,
        EXON_LIST=exons,
        output=output,
        flank=flank,
    )
    gene.run()


def assemble_locus(input: str, blastdb: str, locus: str, output: str, flank: int):
    """
    Assemble a gene from a locus defined in a bed file
    """
    gene = LocusAssembler(
        BED_FILE=input, BLASTDB_PATH=blastdb, locus=locus, output=output, flank=flank
    )
    gene.run()
