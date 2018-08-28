import sys
import argparse
import os
import pandas as pd
from Choruslib.revcom import revcom

def main():

    args = check_options(get_options())


    try:

        probe = pd.read_table(args.input,names=("chrom","start","end","sequence","score","strand"),
                      dtype={"chrom":str,"start":int,"end":int,"sequence":str,"score":int,"strand":str})

        ftprobe = probe[((probe.score > args.mink) & (probe.score < args.maxk) & (
                    probe.score > probe.score.quantile(args.minquantile)) & (
                                 probe.score < probe.score.quantile(args.maxquantile)))]

        probe = ''

        chromnames = ftprobe.chrom.unique()

        outio = open(args.output, 'w')

        for chrom in chromnames:

            prestart = 0

            strandfr = 1

            for index, row in ftprobe[ftprobe.chrom == chrom].sort_values(by='start').iterrows():

                nowpbstart = row.start

                if nowpbstart > prestart + args.dis:

                    if args.strand == True:

                        if strandfr % 2 == 0:
                            row[5] = '-'
                            row[3] = revcom(row[3])

                    strandfr += 1

                    prestart = nowpbstart

                    print(row.chrom, row.start, row.end, row.sequence, row.score, row.strand, sep='\t', file=outio)

        outio.close()

        print("Filtered finished!")

    except:

        print("Probe file format error. Please check your probe file!")

        sys.exit(1)






def check_options(parser):

    args = parser.parse_args()

    if args.input:

        if not os.path.exists(args.input):

            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    if args.mink > args.maxk:

        print("min kmer (%s) score greater than max kmer (%s) score" % (args.mink , args.maxk))

        parser.print_help()

        sys.exit(1)


    if args.minquantile > args.maxquantile:

        print("min quantile (%s) greater than max quantile (%s)" % (args.minquantile , args.maxquantile))

        parser.print_help()

        sys.exit(1)

    if args.minquantile > 1 or args.minquantile < 0:

        print("min quantile error, please check!")

        parser.print_help()

        sys.exit(1)

    if args.maxquantile > 1 or args.maxquantile < 0:

        print("max quantile error, please check!")

        parser.print_help()

        sys.exit(1)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='input file generate by ChorusKmerfilter',
                        required=True, type=str)

    parser.add_argument('-o', '--output', dest='output', help='output file, default=filtered_output.bed',
                        default='filtered_output.bed', type=str)

    parser.add_argument('-m', '--min', dest='mink', help='mini kmer score, example 900, default=0',
                        default=0, type=int)

    parser.add_argument('-l', '--max', dest='maxk', help='max  kmer score, example 2000, default=10000000',
                        default=10000000, type=int)

    parser.add_argument('-q', '--minquantile', dest='minquantile', type=float,
                        help='filter < min\% score range from (0-1)', default=0.1)

    parser.add_argument('-p', '--maxquantile', dest='maxquantile', type=float,
                        help='filter > max\% score range from (0-1)', default=0.9)

    parser.add_argument('-bs', '--bothstrand', dest='strand', help='use both + and - strand, default is True',
                        action='store_true')

    parser.add_argument('-ss', '--singlestrand', dest='strand', help='use only + strand',
                        action='store_false')

    parser.set_defaults(strand=True)

    parser.add_argument('-d', '--distance', dest='dis', help='Min distance between two adjacent probes, default=25',
                        default=25, type=int)

    return parser


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)