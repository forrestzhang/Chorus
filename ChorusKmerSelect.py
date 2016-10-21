import sys
import argparse
import os


def main():

    args = check_options(get_options())

    # input file
    # chrom start   end seq    ChorusKmerfilter_selected   score

    # probes[chrom][start]='chrom\tstart\tend\tseq\tscore'

    probes = dict()

    print(args)

    if args.fast:

        with open(args.input) as inio:

            for lin in inio:

                lintb = lin.rstrip().split('\t')

                if lintb[4] == 'True':

                    if lintb[0] not in probes:

                        probes[lintb[0]] = dict()

                    del lintb[4]

                    nowscore = int(lintb[-1])

                    if args.mink <= nowscore <= args.maxk:

                        probes[lintb[0]][int(lintb[1])] = '\t'.join(lintb)


    else:

        with open(args.input) as inio:

            for lin in inio:

                lintb = lin.rstrip().split('\t')

                if lintb[0] not in probes:

                    probes[lintb[0]] = dict()

                del lintb[4]

                nowscore = int(lintb[-1])

                if args.mink <= nowscore <= args.maxk:

                    probes[lintb[0]][int(lintb[1])] = '\t'.join(lintb)


    outio = open(args.output, 'w')



    for chrom in sorted(probes):

        prestart = 0

        strandfr = 1

        for nowpbstart in sorted(probes[chrom]):

            if nowpbstart > prestart+args.dis:

                nowprob = probes[chrom][nowpbstart].split('\t')

                if args.strand == True:

                    if strandfr % 2 == 0:

                        nowprob.append('-')

                        nowprob[3] = revcom(nowprob[3])

                    else:

                        nowprob.append('+')
                else:

                    nowprob.append('+')

                strandfr += 1

                prestart = nowpbstart

                print('\t'.join(nowprob), file=outio)

    outio.close()

    # out bed format
    # chrom start   end seq    score   strand


def check_options(parser):

    args = parser.parse_args()

    if args.input:

        if not os.path.exists(args.input):

            print("Can not locate input file, please input input file.\n")

            parser.print_help()

            sys.exit(1)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='input file generate by ChorusKmerfilter',
                        required=True, type=str)

    parser.add_argument('-o', '--output', dest='output', help='output file, default=output.bed',
                        default='output.bed', type=str)

    parser.add_argument('-m', '--min', dest='mink', help='mini kmer score, example 900, default=0',
                        default=0, type=int)

    parser.add_argument('-l', '--max', dest='maxk', help='max  kmer score, example 2000, default=10000000',
                        default=10000000, type=int)

    parser.add_argument('-f', '--fast', dest='fast', help='fast filter model, use ChorusKmerSelect filter only',
                        action='store_true')

    parser.add_argument('-nf', '--no-fast', dest='fast', help='fast filter model, DO NOT use ChorusKmerSelect filter',
                        action='store_false')

    parser.set_defaults(fast=True)

    parser.add_argument('-bs', '--bothstrand', dest='strand', help='use both + and - strand, default is True',
                        action='store_true')

    parser.add_argument('-ss', '--singlestrand', dest='strand', help='use only + strand',
                        action='store_false')

    parser.set_defaults(strand=True)

    parser.add_argument('-d', '--distance', dest='dis', help='Min distance between two adjacent probes, default=25',
                        default=25, type=int)

    return parser


def revcom(sequence):

    revdic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','a':'t','t':'a','c':'g','g':'c'}

    return "".join([revdic[base] for base in reversed(sequence)])


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)