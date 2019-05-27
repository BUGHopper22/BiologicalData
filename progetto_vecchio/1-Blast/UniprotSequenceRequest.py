
import os
import sys
import urllib.request

# input: file ids.txt contenente tutti gli id delle proteine
# output: file sequence.fasta con le sequenze delle proteine ()


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
    idList = idProteinList(idFile)
    IdsToSequences(idList)
