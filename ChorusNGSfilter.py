import argparse
import sys
from Choruslib import jellyfish
import os
from multiprocessing import Pool
from pyfasta import Fasta
import pyBigWig

def main():

    args = check_options(get_options())

    # jfgeneratorscount(jfpath, mer, output, generators,threads=1,  size='100M'):

    # make generators

    print(args.input)

    jellyfish.makegenerator(filenames=args.input.split(','), type=args.gzip, generators='generators')

    jfkmerfile = args.output+'.jf'

    bwfile = args.output+'.bw'

    outfilename = args.output

    # jellyfish.jfgeneratorscount(jfpath=args.jellyfish, mer=args.kmer, output=jfkmerfile,
    #                             generators='generators',threads=args.threads,  size='100M')

    spsize = 10000000

    fastain = Fasta(args.genome)

    bw = pyBigWig.open(bwfile, "w")

    seqlenth = dict()
    seqname = dict()

    for chrom in sorted(fastain.keys()):
        infor = chrom.split()
        seqlenth[infor[0]] = len(fastain[chrom])
        seqname[infor[0]] = chrom

    bw.addHeader(list(seqlenth.items()))

    jfscoerlist = list()

    for seqfullname in sorted(fastain.keys()):

        infor = seqfullname.split()

        chrlen = len(fastain[seqfullname])

        if chrlen < spsize:

            start = 0

            end = chrlen - 1

            jfscoer = jellyfish.JFNGSScoer(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=args.kmer,
                                           start=start, end=end, seqfullname=seqfullname, pyfasta=fastain)

            # jfscoer = jellyfish.jfngsscoer(jfscoer)

            #         bw.addEntries(jfscoer.seqname,jfscoer.start,values=jfscoer.score,span=1,step=1)

            jfscoerlist.append(jfscoer)

        else:

            chrblock = int(chrlen / spsize) + 1

            for i in range(chrblock):

                start = i * spsize

                end = start + spsize - 1

                if i > 0:
                    start = start - args.kmer + 1

                if end >= chrlen:
                    end = chrlen - 1

                jfscoer = jellyfish.JFNGSScoer(jfpath=args.jellyfish, jfkmerfile=jfkmerfile, mer=args.kmer,
                                               start=start, end=end, seqfullname=seqfullname, pyfasta=fastain)
                jfscoerlist.append(jfscoer)

    jfsllength = int(len(jfscoerlist)/args.threads + 1)

    for jt in range(jfsllength+1):

        if jt == jfsllength:

            nowlist = jfscoerlist[jt*args.threads+1:]

        else:

            nowlist = jfscoerlist[jt * args.threads:((jt +1) * args.threads -1)]

        pool = Pool(args.threads)

        for jfscoer in pool.map(jellyfish.jfngsscoer, nowlist):

            bw.addEntries(jfscoer.seqname, jfscoer.start, values=jfscoer.score, span=1, step=1)

            print(jfscoer.seqname, jfscoer.start, 'OK')

        pool.close()

        bw.close()


    # pool = Pool(args.threads)
    #
    # for jfscoer in pool.map(jellyfish.jfngsscoer, jfscoerlist):
    #
    #     bw.addEntries(jfscoer.seqname, jfscoer.start, values=jfscoer.score, span=1, step=1)
    #
    #     print(jfscoer.seqname, jfscoer.start, 'OK')
    #
    # bw.close()
    #
    # pool.close()

    bwforcount = pyBigWig.open(bwfile)



    outio = open(outfilename, 'w')

    with open(args.probe) as inio:

        for i in inio:

            (chrom, start, end, seq) = i.rstrip().split()

            score = sum(bwforcount.values(chrom, int(start) - 1, int(end) - 16))

            print(chrom, start, end, seq, int(score), '+', sep='\t', file=outio)

    outio.close()

    print("finished!")

def check_options(parser):

    args = parser.parse_args()

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

    if not os.path.exists(args.genome):

        print("Can not locate genome file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.probe):

        print("Can not locate probe file, please input genome file.\n")

        parser.print_help()

        sys.exit(1)


    if args.input:

        inputfiles = args.input.split(',')

        for inputfile in inputfiles:

            if not os.path.exists(inputfile):

                print("Can not locate %s file.\n" % inputfile)

                parser.print_help()

                sys.exit(1)

    else:

        print("Can not locate input file, please input input file.\n")

        parser.print_help()

        sys.exit(1)

    return args

def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='jellyfish path')

    parser.add_argument('-g', '--genome', dest='genome', help='fasta format genome file', required=True)

    parser.add_argument('-i', '--input', dest='input',
                        help='fastq format input files, spread with \",\", example: 1.fq.gz,2.fq.gz ', required=True,
                        type=str)

    parser.add_argument('-z', '--gziped', dest='gzip', help='gziped file or not, gz or text', default='gz', required=True)

    # parser.add_argument('-s', '--save', dest='saved', help='result saved folder', default='probes')

    parser.add_argument('-t', '--threads', dest='threads', help='threads number or how may cpu you wanna use',
                        default=1, type=int)

    parser.add_argument('-k', '--kmer', dest='kmer', help='kmer length', default=17, type=int)

    parser.add_argument('-p', '--probe', dest='probe', help='probe file')

    parser.add_argument('-o', '--output', dest='output', help='output file', default='output.bed')

    # args = parser.parse_args()

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
