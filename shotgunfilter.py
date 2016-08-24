import os
import subprocess
import pandas as pd
import argparse
import sys
from Choruslib import jellyfish


def main():

    args = check_options(get_options())

    pbfile = args.input

    outfile = args.outfile

    jfkmerfile = args.kmerfile

    maxcount = args.maxcount

    mincount = args.mincount

    jfpath = args.jellyfish

    mer = args.mer

    pbdf = pd.read_table(pbfile,names=("Chr","start","end","seq"))

    sequencelist = pbdf.seq.tolist()

    ct = jfseqkmertestshotgun(jfpath, jfkmerfile, mer, sequencelist, maxcount, mincount)

    pbdf[ct].to_csv(outfile, sep='\t', header=False, index=False)

#
#     pbfile = 'Zea_mays.AGPv3.26.dna_rm.genome.fa_all.bed'
#
#     outfile = 'shotgunfilter_'+pbfile
#
#     pbdf = pd.read_table(pbfile,names=("Chr","start","end","seq"))
#
#     test = pbdf.seq.tolist()
#
#     ct = jfseqkmertestshotgun('/opt/bin/jellyfish', '/home/forrest/FISH/maize/SRA/test_maize17mer.jf', 17, test, 320)
#
#     pbdf[ct].to_csv(outfile,sep='\t',header=False,index=False)


def subprocesspath(path):
    """

    :param path: path
    :return: path, for subprocess, avoid white space error
    """

    rpath = '\''+ os.path.abspath(path)+'\''

    return rpath


def jfseqkmertestshotgun(jfpath, jfkmerfile, mer, sequencelist, maxcount, mincount):

    """
    :param jfpath: jellyfish bin path
    :param jfkmerfile: jellyfish kmer count file
    :param mer: int, kmer
    :param sequence: string, sequence for kmerscore count
    :param bfcount:
    :return: list, kmerscore list
    """

    jfpath = subprocesspath(jfpath)

    jfkmerfile = subprocesspath(jfkmerfile)

    jfquerycommand = ' '.join([jfpath, 'query', '-i', '-l', jfkmerfile])

    print(jfquerycommand)

    kmerct = subprocess.Popen(jfquerycommand, shell=True, stdout=subprocess.PIPE,
                                  stdin=subprocess.PIPE)

    filterlist = list()

    for sequence in sequencelist:
        mer = int(mer)

        end = mer
        seqlen = len(sequence)

        kmerfilter = True

        while (end <= seqlen):

            start = end - mer

            subseq = sequence[start:end]+'\n'

            kmerct.stdin.write(subseq.encode('ascii'))

            kmerct.stdin.flush()

            lin = kmerct.stdout.readline().decode('utf-8').rstrip('\n')

            number = int(lin)

            if number > maxcount:

                kmerfilter = False

            if number < mincount:

                kmerfilter = False

            end +=1

        filterlist.append(kmerfilter)

    kmerct.stdin.close()

    kmerct.stdout.close()

    kmerct.wait()

    return filterlist



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

    if not os.path.exists(args.input):

        print("Can not input file, please input probe file.\n")

        parser.print_help()

        sys.exit(1)

    if args.mincount > args.maxcount:

        print("min count %s is larger than max count %s" % (args.mincount, args.maxcount))

        parser.print_help()

        sys.exit(1)

    print("jellyfish version:", args.jellyfish, jellyfishversion)

    print("input file:", args.input)

    print("min count %s, max count %s" % (args.mincount, args.maxcount))

    print("result output folder:", os.path.realpath(args.outfile))

    return args


def get_options():

    parser = argparse.ArgumentParser(description="Chorus Software for Oligo FISH probe design", prog='Chorus')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-i', '--input', dest='input', help='unfiltered prob file, generate by Chrus.', required=True)

    parser.add_argument('-k', '--kmer', dest='kmerfile', help='kmerfile', required=True)

    parser.add_argument('--mer', dest='mer', default=17, type=int)

    parser.add_argument('-j', '--jellyfish', dest='jellyfish', help='jellyfish path')

    parser.add_argument('-o', '--out', dest='outfile', help='output file name', default='outprobe.txt')

    parser.add_argument('-l', '--min', dest='mincount', help='min count of kmer', type=int, default=2)

    parser.add_argument('-m', '--max', dest='maxcount', help='max count of kmer', type=int, required=True)

    return parser

def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates





if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


            #
# if __name__ == "__main__":
#
#     pbfile = 'Zea_mays.AGPv3.26.dna_rm.genome.fa_all.bed'
#
#     outfile = 'shotgunfilter_'+pbfile
#
#     pbdf = pd.read_table(pbfile,names=("Chr","start","end","seq"))
#
#     test = pbdf.seq.tolist()
#
#     ct = jfseqkmertestshotgun('/opt/bin/jellyfish', '/home/forrest/FISH/maize/SRA/test_maize17mer.jf', 17, test, 320)
#
#     pbdf[ct].to_csv(outfile,sep='\t',header=False,index=False)