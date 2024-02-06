import pandas as pd
import click
import sys
import logging


class BlastConverter:
    def __init__(
        self,
        input: str,
        locus_size: int,
        exon_count: int,
        q_cov_threshold: float,
        refseq=False,
    ) -> None:
        self.input_file = input
        self.cds_blast_data = pd.read_csv(
            self.input_file,
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
        self.bed9 = pd.DataFrame()
        self.locus_size = locus_size
        self.exon_count = exon_count
        self.q_cov_threshold = q_cov_threshold

        logging.info(f"Converting {self.input_file} to bed format...")
        logging.info(f"# of Exons: {self.exon_count}")
        logging.info(f"Locus size: {self.locus_size}")
        logging.info(f"Coverage threshold: {self.q_cov_threshold}")

    def process_BLAST(self):
        # Establish strand orientation and query coverage of BLAST hits
        self.cds_blast_data["orientation"] = (
            self.cds_blast_data.sstart < self.cds_blast_data.send
        )
        self.cds_blast_data["strand"] = self.cds_blast_data.orientation.map(
            lambda x: "+" if x is True else "-"
        )
        self.cds_blast_data["qcov"] = round(
            self.cds_blast_data.length / self.cds_blast_data.qlen, 2
        )

        # Correctly order start/end of BLAST hits for BED
        cond = self.cds_blast_data.sstart > self.cds_blast_data.send
        self.cds_blast_data.loc[cond, ["sstart", "send"]] = self.cds_blast_data.loc[
            cond, ["send", "sstart"]
        ].values

        # Sort values + reindex
        self.cds_blast_data.sort_values(by="sstart", inplace=True)
        self.cds_blast_data.reset_index(drop=True, inplace=True)

        self.cds_blast_data["chromosome"]: str = self.cds_blast_data.sseqid.map(lambda x: x.split("|")[1] if "|" in x else x) 


    def convert_BLAST_to_BED(self):
        # Fill columns 1-9
        self.bed9["chrom"]: str = self.cds_blast_data.chromosome  # type: ignore
        self.bed9["chromStart"]: int = self.cds_blast_data.sstart  # type: ignore
        self.bed9["chromEnd"]: int = self.cds_blast_data.send  # type: ignore
        self.bed9["name"]: str = [f"exon_{i}" for i in self.cds_blast_data.index]  # type: ignore
        self.bed9["score"]: float = self.cds_blast_data.qcov  # type: ignore
        self.bed9["strand"]: str = self.cds_blast_data.strand  # type: ignore
        self.bed9["thickStart"]: int = self.cds_blast_data.sstart  # type: ignore
        self.bed9["thickEnd"]: int = self.cds_blast_data.send  # type: ignore
        self.bed9["itemRgb"]: str = "145,30,180"  # type: ignore
        self.bed9["exonList"]: str = ""

    def predict_gene_loci(self):
        gene_locus = 0

        for window in self.cds_blast_data.sort_values(["chromosome", "strand"]).rolling(
            self.exon_count
        ):
            if (
                len(window) > self.exon_count - 1
                and all(window.chromosome.unique())  # type: ignore
                and all(window.strand.unique())
            ):
                locus_qcov = round(window.qcov.sum(), 2)
                size = window.send.max() - window.sstart.min()

                if locus_qcov > self.q_cov_threshold and size < self.locus_size * 1.5:
                    self.bed9.loc[len(self.bed9.index)] = [
                        window.chromosome.unique()[0],
                        window.sstart.min(),
                        window.send.max(),
                        f"locus_{gene_locus}",
                        locus_qcov,
                        window.strand.unique()[0],
                        window.sstart.min(),
                        window.send.max(),
                        "0,255,0",
                        ""
                    ]
                    gene_locus += 1
        logging.info(f"{gene_locus} gene loci predicted...")

    def write_bed_file(self):
        # Save output to bed file (formatted as tsv)
        self.bed9.to_csv(
            sys.stdout,
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

    def run(self):
        self.process_BLAST()
        self.convert_BLAST_to_BED()
        self.predict_gene_loci()
        self.write_bed_file()


@click.command()
@click.option(
    "-i",
    "--input",
    type=click.Path(exists=True),
    required=True,
    help="BLAST file in tabular format",
)
@click.option(
    "-e", "--exons", type=int, default=0, help="Expected number of exons in the gene"
)
@click.option(
    "-c",
    "--coverage",
    type=float,
    default=1.1,
    help="Proportion of gene that must be covered by a predicted locus",
)
@click.option(
    "-l", "--locus_size", type=int, default=1000, help="Expected size of the locus"
)
def blast2bed(input: str, exons: int, coverage: int, locus_size: int):
    """
    Convert BLAST results to BED and predict gene loci
    """
    converter = BlastConverter(
        input=input,
        exon_count=exons,
        q_cov_threshold=coverage,
        locus_size=locus_size
    )
    converter.run()
