import click
from blast_2_bed import blast2bed
from assemble_gene import GeneAssembler

@click.group(
    "S.A.B.A.T: Semi-Automatic BLAST Annotation Toolkit"
)
def cli() -> None:
    pass

cli.add_command(blast2bed)

if __name__ == "__main__":
    cli()