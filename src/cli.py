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


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")
    
    else:
        return arg


def main():
    # Initialise parser and subparser for each command
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="command", prog="sabat")

    # Parser for blast2bed command 
    blast2bed_parser = sub_parsers.add_parser("blast2bed", help="Convert BLAST output to BED and predict gene loci")

    # Blast2bed arguments
    blast2bed_parser.add_argument('input', metavar="input", type=lambda x: is_valid_file(parser, x), help="BLAST file in tabular format (required)")
    blast2bed_parser.add_argument('-e', "--exons", metavar="exons", type=int, required=False, default=0, help="Expected number of exons in the gene (default: 0)")
    blast2bed_parser.add_argument("-c", "--coverage", metavar="coverage", type=float, required=False, default=1.1, help="Proportion of gene that must be covered by a predicted locus (default: 1.1)")
    blast2bed_parser.add_argument("-l", "--locus_size", metavar="locus_size", type=int,required=False, default=1000, help="Expected size of the locus in base pairs (default: 1000)")

    # Assemble-exons parser
    assemble_exons_parser = sub_parsers.add_parser("assemble-exons", help="Assemble gene from specified exons")
    
    # Assemble-exons arguments
    assemble_exons_parser.add_argument("-i", "--input", metavar="input", type=lambda x: is_valid_file(parser, x), required=True, help="bed file")
    assemble_exons_parser.add_argument("-db", "--blastdb", metavar="blastdb", required=True, help="Path to the BLASTdb (including name)")
    assemble_exons_parser.add_argument("-e", "--exons", metavar="exons", nargs='+', type=str, required=True, help="Name of exons")
    assemble_exons_parser.add_argument("-o", "--output", metavar="output", required=True, help="Base name for output files")
    assemble_exons_parser.add_argument("-f","--flank", metavar="flank", type=int, default=0, required=False, help="Number of nucleotides to add to 3' flank of the predicted gene (ensures stop codon is found)")

    # Assemble-locus parser
    assemble_locus_parser = sub_parsers.add_parser("assemble-locus", help="Assemble gene from a defined locus")

    # Assemble-locus arguments
    assemble_locus_parser.add_argument("-i", "--input", metavar="input", type=lambda x: is_valid_file(parser, x), required=True, help="bed file")
    assemble_locus_parser.add_argument("-db", "--blastdb", metavar="blastdb", required=True, help="Path to the BLASTdb (including name)")
    assemble_locus_parser.add_argument("-l", "--locus", metavar="locus", type=str, required=True, help="Name of the predicted locus")
    assemble_locus_parser.add_argument("-o", "--output", metavar="output", required=True, help="Base name for output files")
    assemble_locus_parser.add_argument("-f","--flank", metavar="flank", type=int, default=0, required=False, help="Number of nucleotides to add to 3' flank of the predicted gene (ensures stop codon is found)")

    
    args = parser.parse_args()
    
    if args.command == "blast2bed":
        blast2bed(input=args.input, exons=args.exons, coverage=args.coverage, locus_size=args.locus_size)

    if args.command == "assemble-locus":
        assemble_locus(input=args.input, blastdb=args.blastdb, locus=args.locus, output=args.output, flank=args.flank)
    
    if args.command == "assemble-exons":
        assemble_exons(input=args.input, blastdb=args.blastdb, exons=args.exons, output=args.output, flank=args.flank)




if __name__ == "__main__":
    main()

