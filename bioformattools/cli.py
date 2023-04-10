import argparse
from bioformattools.pipelines import *


def args_init():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='BioFormatTools',
    )

    subparsers = parser.add_subparsers(
        title='subcommands', dest="subcommand_name")

    # argparse for Gff2Gtf
    parser_a = subparsers.add_parser('Gff2Gtf',
                                     help='convert gff to gtf file')

    parser_a.add_argument('input_file', type=str,
                          help='input file with gff format')
    parser_a.add_argument('output_file', type=str, help='output file path')

    # argparse for Gtf2Gff
    parser_a = subparsers.add_parser('Gtf2Gff',
                                        help='convert gtf to gff file')

    parser_a.add_argument('input_file', type=str,
                          help='input file with gtf format')
    parser_a.add_argument('output_file', type=str, help='output file path')

    # argparse for genblasta2BED
    parser_a = subparsers.add_parser('genblasta2BED',
                                     help='convert genblasta output to bed file')

    parser_a.add_argument('input_file', type=str,
                          help='input file with outfmt 5')
    parser_a.add_argument('output_file', type=str, help='output file path')
    parser_a.add_argument('-p', "--ID_prefix", type=str, default='subject_',
                          help='gene output name prefix defaults: subject_')

    # argparse for outfmt5To6
    parser_a = subparsers.add_parser('outfmt5To6',
                                     help='convert blast results from outfmt 5 to 6',
                                     description='convert blast results from outfmt 5 to 6')

    parser_a.add_argument('input_file', type=str,
                          help='input file with outfmt 5')
    parser_a.add_argument('output_file', type=str, help='output file path')

    # argparse for outfmt5complete
    parser_a = subparsers.add_parser('outfmt5complete',
                                     help='check if outfmt5 is complete')

    parser_a.add_argument('input_file', type=str,
                          help='input file with outfmt 5')

    # argparse for outfmt6ToFasta
    parser_a = subparsers.add_parser('outfmt6ToFasta',
                                     help='extract subject sequence by blast outfmt6 results')

    parser_a.add_argument('outfmt6', type=str, help='input file with outfmt 6')
    parser_a.add_argument('db_fasta', type=str,
                          help='input file with database fasta file')
    parser_a.add_argument('output_fasta', type=str, help='output file path')

    args = parser.parse_args()

    return args


def main():

    args = args_init()
    args_dict = vars(args)

    if args_dict["subcommand_name"] == "Gff2Gtf":
        Gff2Gtf(args.input_file, args.output_file)
    elif args_dict["subcommand_name"] == "Gtf2Gff":
        Gtf2Gff(args.input_file, args.output_file)
    elif args_dict["subcommand_name"] == "genblasta2BED":
        genblasta2BED_main(args)
    elif args_dict["subcommand_name"] == "outfmt5To6":
        bf = BlastFormat(input_blast_file=args.input_file, output_file=args.output_file)
        bf.outfmt5To6()
    elif args_dict["subcommand_name"] == "outfmt5complete":
        bf = BlastFormat(input_blast_file=args.input_file)
        bf.outfmt5complete()
    elif args_dict["subcommand_name"] == "outfmt6ToFasta":
        bf = BlastFormat(input_blast_file=args.outfmt6, input_fasta_file=args.db_fasta, output_file=args.output_fasta)
        bf.outfmt6ToFasta()

if __name__ == '__main__':
    main()