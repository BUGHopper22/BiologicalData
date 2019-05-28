
# input: R03_msa200.clw contains the multiple sequence allignement of all sequences
# output: R04_msa200Cleaned.clw
# function: remove row and column with too many gaps

from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
from datetime import date
import os
#import matplotlib.pyplot as plt


def cut_alignment(alignment, threshold_c=0, threshold_r=0):
    """
        the variable threshold_c is the maximum number of dashes accepted in the columns
        the variable threshold_r is the maximum number of dashes accepted in the rows
    """
    l = np.array([np.array(list(str(s.seq))) for s in alignment])

    c = np.sum(l == "-", axis=0)/l.shape[0]

    if threshold_c == 0:
        threshold_c = 0.1

    idx = np.arange(c.shape[0])
    idx = idx[c < threshold_c]

    #c_x = c[idx]
    # plt.plot(idx,c_x)
    new_l = l[:, idx]

    r = np.sum(new_l == "-", axis=1)/new_l.shape[1]

    if threshold_r == 0:
        threshold_r = 0.06

    idx = np.arange(r.shape[0])
    idx = idx[r > threshold_r]

    new_l[idx, :] = '?'  # check to not include undesired sequences

    seq_n = [None]*len(alignment)

    for i in range(len(alignment)):
        if i not in idx:
            s = list(new_l[i, :])
            s = [str(el) for el in s]
            seq_n[i] = "".join(s)

    msa_tmp = [None]*len(alignment)
    msa = []

    for i in range(len(alignment)):
        if i not in idx:
            msa_tmp[i] = SeqRecord(seq=Seq(seq_n[i], SingleLetterAlphabet()), id=alignment[i].id,
                                   name=alignment[i].name, description=alignment[i].description, dbxrefs=alignment[i].dbxrefs)
            msa.append(msa_tmp[i])

    align = MultipleSeqAlignment(msa)
    return align


# __CHECK IF ALREADY EXIST FILENAME
# _return true iif already exist filename
def alreadyExistNameFile(folderNameToSave, path):
    # filenameToSave become unique
    for root, dirs, files in os.walk(path):
        trovato = False
        for folder in dirs:
            folder = folder + "/"
            if(folder == folderNameToSave):
                return True
    return False


def FindOutputPathFolder(path):
    currentAttempt = 1
    titleAttempt = '{0:03}'.format(currentAttempt)
    folderNameToSave = str(currentDate) + "_" + str(titleAttempt) + "/"

    # if name exist increment titleAttempt
    while(alreadyExistNameFile(folderNameToSave, path)):
        currentAttempt += 1
        currentTitleAttempt = '{0:03}'.format(currentAttempt)
        folderNameToSave = str(currentDate) + "_" + \
            str(currentTitleAttempt) + "/"

    outputFolderPathName = path + folderNameToSave

    return outputFolderPathName


if __name__ == "__main__":

    currentDate = today = date.today()
    d1 = today.strftime("%Y-%m-%d")

    # __INPUT
    inputFilePath = "../results/R03_msa200.clw"
    fileToAlign = AlignIO.read(open(inputFilePath), "clustal")

    # ____________ CHANGE THE RANGES HERE! ____________

    beginThresholdCol = 0.20
    endThresholdCol = 0.60

    beginThresholdRow = 0.20
    endThresholdRow = 0.70

    # ____________ CHANGE THE FOLDERPATH HERE! ____________
    path = "../results/R04_msaCleanedTest/"
    # _________________________________________________
    outputPathFolder = FindOutputPathFolder(path)
    bufferingCount = 1
    for actualRow in np.arange(beginThresholdRow, endThresholdRow, 0.05):

        for actualCol in np.arange(beginThresholdCol, endThresholdCol, 0.05):

            alignmentOutput = cut_alignment(fileToAlign, actualCol, actualRow)

            if not os.path.exists(outputPathFolder):
                os.makedirs(outputPathFolder)

            # title with threshold
            titleCol = actualCol*100
            titleRow = actualRow*100

            outputFileNamePath = outputPathFolder + "msa_row" + "{:.0f}".format(titleCol) + \
                "_col" + "{:.0f}".format(titleRow) + ".clw"

            # stampa di attesa:
            numberOfFiles = (beginThresholdCol + endThresholdCol)*20
            bufferingCount += 1
            print("81 files, {} completed".format(bufferingCount))

            outputFile = open(outputFileNamePath, "w")

            AlignIO.write(alignmentOutput, outputFile, "clustal")
            outputFile.close()
    print("FINISH")
