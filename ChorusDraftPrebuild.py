from pyfasta import Fasta
import sys
import argparse
import os


def check_options(parser):

    args = parser.parse_args()

    if args.input:

        if not os.path.exists(args.input):

            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="Combine short sequence to speed oligo search")

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='input file',
                        required=True, type=str)

    parser.add_argument('-o', '--output', dest='output', help='output file, default=output.fa',
                        default='output.fa', type=str)


def main():

    args = check_options(get_options())

    fain = Fasta(args.input)

    faout = open(args.output, 'w')

    minlen = 1e6

    shortseq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

    breacker = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

    for chrome in fain.keys():

        if len(fain[chrome]) < minlen:
            shortseq = shortseq + fain[chrome] + breacker

        else:

            print('>%s' % chrome, file=faout)
            print(fain[chrome], file=faout)

    print('>shortsequences', file=faout)

    print(shortseq, file=faout)


    faout.close()


if __name__ == "main":

    main()
