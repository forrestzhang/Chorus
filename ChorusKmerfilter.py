from Choruslib import jellyfish
import sys
import argparse
import os
from Chorus import which
from multiprocessing import Pool


def main():

    args = check_options(get_options())

    jfpool = Pool(args.threads)

    probelist = list()

    with open(args.input, 'r') as probein:
        for i in probein:

            jfpbkfruner = jellyfish.JFpbkfruner(jfpath=args.jellyfish, jfkmerfile=args.kmerfile,mink=args.mink, maxk=args.maxk, mer=args.kmer, probe=i.strip('\n'))

            probelist.append(jfpbkfruner)

    probein.close()

    dirout = os.path.dirname(args.input)

    basename = os.path.basename(args.input)

    outfilename = os.path.join( dirout,'kmerfilted_'+basename)

    probeout = open(outfilename, 'w')

    for currentprobe in jfpool.imap_unordered(jellyfish.jfprobekmerfilter, probelist):

        print(currentprobe['chro'] ,
                currentprobe['start'] ,
                currentprobe['end'] ,
                currentprobe['seq'] ,
                currentprobe['keep'] ,
                currentprobe['sumscore'],file=probeout,sep='\t')

    probeout.close()



def check_options(parser):

    args = parser.parse_args()

    # Start check jellyfish
    if args.jellyfish:

        if not os.path.exists(args.jellyfish):

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)

        jellyfishversion = jellyfish.jfversion(args.jellyfish)

        if jellyfishversion == 'None':

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)

    else:

        jellyfishpath = which('jellyfish')

        if jellyfishpath:

            jellyfishversion = jellyfish.jfversion(jellyfishpath[0])

            if jellyfishversion == 'None':

                print("Can not locate jellyfish, please input full path of jellyfish\n")

                parser.print_help()

                sys.exit(1)

            else:

                args.jellyfish = jellyfishpath[0]

        else:

            print("Can not locate jellyfish, please input full path of jellyfish\n")

            parser.print_help()

            sys.exit(1)
    # End check jellyfish


    # check max min
    if args.maxk < args.mink:

        print("Error! Max kmer count is smaller than min kmer count!\n")

        parser.print_help()

        sys.exit(1)
    # End check max min

    # Check input file
    if not os.path.exists(args.input):

        print("Can not locate input file, please input input file.\n")

        parser.print_help()

        sys.exit(1)
    # End check input file
    return args


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='input file generate by Chorus',
                        required=True, type=str)

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='jellyfish path')

    parser.add_argument('-t', '--threads', dest='threads', help='threads number or how may cpu you wanna use',
                        default=1, type=int)

    parser.add_argument('-k', '--kmer', dest='kmer', help='kmer length',
                        default=17, type=int)

    parser.add_argument('-m', '--min', dest='mink', help='mini kmer count, example 10',
                        required=True, type=int)

    parser.add_argument('-l', '--max', dest='maxk', help='max kmer count, example 50',
                        required=True, type=int)

    parser.add_argument('-f', '--kmerfile', dest='kmerfile', help='kmer file',
                        required=True, type=str)

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)