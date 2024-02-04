import logging
import sys
import click
from blast_2_bed import blast2bed
from assemble_gene import assemble_gene

logging.basicConfig(
    stream=sys.stderr,
    format="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.INFO,
)


@click.group(help="S.A.B.A.T: Semi-Automatic BLAST Annotation Toolkit")
def cli() -> None:
    pass


cli.add_command(blast2bed)
cli.add_command(assemble_gene)

if __name__ == "__main__":
    cli()
