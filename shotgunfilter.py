__author__ = 'Forrest'

import os
import subprocess
import pandas as pd

def subprocesspath(path):
    """

    :param path: path
    :return: path, for subprocess, avoid white space error
    """

    rpath = '\''+ os.path.abspath(path)+'\''

    return rpath


def jfseqkmertestshotgun(jfpath, jfkmerfile, mer, sequencelist, threshold, bfcount=False):

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




    jfkmercount = list()

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

            if number > threshold:

                kmerfilter = False

#             print(start)

            end +=1

        filterlist.append(kmerfilter)

    kmerct.stdin.close()

    kmerct.stdout.close()

        # kmerct.terminate()

    kmerct.wait()

    return filterlist


if __name__ == "__main__":

    pbfile = 'Zea_mays.AGPv3.26.dna_rm.genome.fa_all.bed'

    outfile = 'shotgunfilter_'+pbfile

    pbdf = pd.read_table(pbfile,names=("Chr","start","end","seq"))

    test = pbdf.seq.tolist()

    ct = jfseqkmertestshotgun('/opt/bin/jellyfish', '/home/forrest/FISH/maize/SRA/test_maize17mer.jf', 17, test, 320)

    pbdf[ct].to_csv(outfile,sep='\t',header=False,index=False)