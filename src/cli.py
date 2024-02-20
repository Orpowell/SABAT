import logging
import sys
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

def main():
    parser = argparse.ArgumentParser()

    sub_parsers = parser.add_subparsers(dest="command")

    blast2bed_parser = sub_parsers.add_parser("blast2bed", help="Convert BLAST output to BED and predict gene loci")


    exon_parser = sub_parsers.add_parser("assemble-exons", help="Assemble gene from specified exons")
    exon_parser.add_argument("-i", metavar="--input", help="path to bed file")

    locus_parser = sub_parsers.add_parser("assemble-locus", help="Assemble gene from a defined locus")

    args = parser.parse_args()

    print(args)




if __name__ == "__main__":
    main()
    #cli()
