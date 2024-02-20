import logging
import sys
import os
import argparse
from blast_2_bed import blast2bed
from assemble_gene import assemble_exons
from assemble_gene import assemble_locus

logging.basicConfig(
    stream=sys.stderr,
    format="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.INFO,
)
'''
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
'''
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")
    
    else:
        return arg


def main():
    parser = argparse.ArgumentParser()

    sub_parsers = parser.add_subparsers(dest="command", prog="sabat")

    # Parser for blast2bed command 
    blast2bed_parser = sub_parsers.add_parser("blast2bed", help="Convert BLAST output to BED and predict gene loci")

    # Blast2bed arguments
    blast2bed_parser.add_argument('-i', "--input", metavar="input", type=lambda x: is_valid_file(parser, x), required=True, help="BLAST file in tabular format (required)")
    blast2bed_parser.add_argument('-e', "--exons", metavar="exons", type=int, required=False, default=0, help="Expected number of exons in the gene (default: 0)")
    blast2bed_parser.add_argument("-c", "--coverage", metavar="coverage", type=float, required=False, default=1.1, help="Proportion of gene that must be covered by a predicted locus (default: 1.1)")
    blast2bed_parser.add_argument("-l", "--locus_size", metavar="locus_size", type=int,required=False, default=1000, help="Expected size of the locus in base pairs (default: 1000)")

    # Assemble-exons parser
    assemble_exons_parser = sub_parsers.add_parser("assemble-exons", help="Assemble gene from specified exons")
    
    # Assemble-exons arguments
    assemble_exons_parser.add_argument("-i", metavar="--input", help="path to bed file")

    # Assemble-locus parser
    assemble_locus_parser = sub_parsers.add_parser("assemble-locus", help="Assemble gene from a defined locus")

    # Assemble-locus arguments
    #TBA
    
    args = parser.parse_args()
    
    if args.command == "blast2bed":
        blast2bed(input=args.input, exons=args.exons, coverage=args.coverage, locus_size=args.locus_size)

    




if __name__ == "__main__":
    main()
    #cli()
