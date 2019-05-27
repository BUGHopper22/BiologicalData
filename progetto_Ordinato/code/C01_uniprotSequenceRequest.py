
# input: R01_idsFromBlast.txt contains all ids of proteins (NCBI BLAST+)
# output: R02_sequence.fasta contains all the sequences of ids
# function: search in uniprot all the sequences of ids.

import os
import sys
import urllib.request


def idProteinList(idFile):
    return [l[15:-1] for l in open(idFile)]


def IdsToSequences(idList):
    uniprotUrl = 'https://www.uniprot.org/uniprot/'
    upiparcUrl = 'https://www.uniprot.org/uniparc/'
    # extension = '.fasta'
    with open('sequences.fasta', 'w') as output:
        for i in range(0, 200):
            if(idList[i][:3] != 'UPI'):
                url = uniprotUrl + idList[i] + '.fasta'
                singleProtein = urllib.request.urlopen(url).read()
                output.write(singleProtein.decode("utf-8"))
            else:
                url = upiparcUrl + idList[i] + '.fasta'
                singleProtein = urllib.request.urlopen(url).read()
                output.write(singleProtein.decode("utf-8"))


if __name__ == "__main__":
    idFile = 'ids.txt'
    # attenzione funziona solamente su questo DB che abbiamo scelto
    idList = idProteinList(idFile)
    IdsToSequences(idList)
