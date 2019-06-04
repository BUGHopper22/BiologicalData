# %%
from Bio import SearchIO
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import SeqIO
import pandas as pd
import os
import _pickle as pickle
import math


def insertDataInArrays(qresults, scores, group, family):
    for i, qresult in enumerate(qresults.hits):
        scores[i] = qresult.bitscore
        desc = qresult.description.split()
        group[i] = desc[1]
        family[i] = desc[2]

# create three arrays with statistical measures for each possible threshold


def buildPerformanceMeasuresLists(scores, myFamilyTrue):

    # last x-score in plot
    maxBitScore = int(np.max(scores)+10)

    # step for score evaluation
    minStep = np.min(-np.diff(scores)) - 0.00001
    if minStep <= 0:
        minStep = 0.09

    # DA CAMBIARE ASSOLUTAMENTE
    # l = int(maxBitScore/minStep)
    l = len(np.arange(0, maxBitScore, minStep))

    accuracies = [0]*l
    sensitivities = [0]*l
    specificities = [0]*l

    bestAcc = -1
    bestAccIndex = -1

    for i, j in enumerate(np.arange(0, maxBitScore, minStep)):

        myFamilyPred = np.int32(scores > j)

        cm = confusion_matrix(myFamilyTrue, myFamilyPred)

        tn = cm[0, 0]
        fp = cm[0, 1]
        fn = cm[1, 0]
        tp = cm[1, 1]

        acc = (tp + tn) / (tp + fp + fn + tn)
        sens = tp / (tp + fn)
        spec = tn / (tn + fp)

        accuracies[i] = acc
        sensitivities[i] = sens
        specificities[i] = spec

        if bestAcc < acc:
            bestAcc = acc
            bestAccIndex = myFamilyPred.sum()-1
            # (tp + tn)
            # len(myFamilyPred == 1) - 1
    return accuracies, sensitivities, specificities, bestAccIndex, bestAcc


def insertDataInKin_ids(kin_ids, thresholdIndex, qresults):
    i = 0
    kin_ids = []

    # insert sequence ids in kin_ids iif bitscore > thresholdIndex
    while qresults.hits[i].bitscore > thresholdIndex:
        # print('qresults.hits[i].bitscore')
        # print(qresults.hits[i].bitscore)
        desc = qresults.hits[i].description
        kin_id = desc.split()[0]
        kin_ids.append(kin_id)
        i = i + 1
    # print(kin_ids)
    return kin_ids


def initializeAllValueTitle(family, myGroupIndexies, familyUniqueList):

    # print('familyUniqueList')
    # print(familyUniqueList)

    allValueTitle = []
    allValueTitle.append("fileName")
    allValueTitle.append("goodnessValue")
    allValueTitle.append("thresholdIndex")
    allValueTitle.append("lenKinIds")
    allValueTitle.append("sensitivityValue")

    # create a List with columns title for initialize resultMatrix
    for i in familyUniqueList:
        allValueTitle.append(i)
    allValueTitle = np.array(allValueTitle)

    return allValueTitle


def analyse_hmm_out(out_hmm, seq_rec_file, resultsMatrix):

    # qResults is a matrix with ID - HSP - DESCRIPTION of every protein in file.tblout
    qresults = next(SearchIO.parse(out_hmm, 'hmmer3-tab'))

    # initialize arrays: scores, group, family
    N = len(qresults.hits)
    scores = [0]*N
    group = [None]*N
    family = [None]*N
    insertDataInArrays(qresults, scores, group, family)

    # convert lists in array
    group = np.array(group)
    scores = np.array(scores)
    family = np.array(family)

    # ______________________CHANGE FAMILY HERE___________________________
    myGroup = 'STE'
    # ___________________________________________________________________

    # boolean matrix where myGroup becomes true
    myGroupBoolArray = (group == myGroup)

    # array with indexies of myGroup
    myGroupIndexies = np.arange(N)[myGroupBoolArray]

    # true become 1, false become 0
    myGroupTrue = np.int32(myGroupBoolArray)

    # create three list with statistical value needed
    accuracies, sensitivities, specificities, bestAccIndex, bestAcc = buildPerformanceMeasuresLists(
        scores, myGroupTrue)
    # print(accuracies)

    thresholdIndex = bestAccIndex
    # print('thresholdIndex')
    print("PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP")
    print(thresholdIndex)
    print(bestAcc)
    print(len(accuracies))
    print(np.argmax(accuracies))
    print(max(accuracies))

    kin_ids = []
    kin_ids = insertDataInKin_ids(kin_ids, thresholdIndex, qresults)

    fileName = out_hmm
    goodnessValue = np.sum((myGroupIndexies + 1)**2)

    lenKinIds = len(kin_ids)
    sensitiveValue = sensitivities[thresholdIndex]
    # print('sensitiveValue')
    # print(sensitiveValue)

    # list of all the family of STE
    familyUniqueList = []
    familyUniqueList = np.unique(family[myGroupIndexies])
    familyUniqueList = np.array(familyUniqueList)
    allValueTitle = initializeAllValueTitle(
        family, myGroupIndexies, familyUniqueList)

    # initialize result matrix dinamically with all title in allValueTitle
    # EQUIVALENT TO resultsMatrix = {"fileName": [], "goodnessValue": [], "thresholdIndex": [
    # ], "lenKinIds": [], "sensitivityValue": [], "valSTE20": [], "valSTE7": [], "valSTE20": []}
    # count = 0
    # resultsMatrix = {}
    # for title in allValueTitle:
    #     resultsMatrix[title] = []
    # print(allValueTitle)
    print(family)
    print(familyUniqueList)
    familyRangeIndexMatrix = {"familyName": [], 'BI': [], 'WI': [], 'SV': []}
    for fam in familyUniqueList:
        familyRangeIndexMatrix["familyName"].append(fam)
        familyIndexTemp = np.arange(N)[family == fam]
        # index of first element with family==fam
        familyRangeIndexMatrix["BI"].append(familyIndexTemp[0])
        # index of last element with family==fam
        familyRangeIndexMatrix["WI"].append(familyIndexTemp[-1])
        totActualFamily = len(familyIndexTemp)
        captureActualFamily = 0
        for index in familyIndexTemp:
            if(index <= thresholdIndex):
                captureActualFamily += 1
        familyRangeIndexMatrix["SV"].append(
            captureActualFamily/totActualFamily)
    print(familyRangeIndexMatrix)

    familyRangeIndexMatrix = {}
    for fam in familyUniqueList:

        totActualFamily = len(familyIndexTemp)
        captureActualFamily = 0
        for index in familyIndexTemp:
            if(index <= thresholdIndex):
                captureActualFamily += 1

        familyRangeIndexMatrix[fam] = []
        familyRangeIndexMatrix[fam].append(captureActualFamily/totActualFamily)

        familyIndexTemp = np.arange(N)[family == fam]
        # index of first element with family==fam
        familyRangeIndexMatrix[fam + "-BI"] = []
        familyRangeIndexMatrix[fam + "-BI"].append(familyIndexTemp[0])
        # index of last element with family==fam
        familyRangeIndexMatrix[fam + "-WI"] = []
        familyRangeIndexMatrix[fam + "-WI"].append(familyIndexTemp[-1])

    fileName = os.path.basename(out_hmm)
    resultsMatrix = {"fileName": [], "goodnessValue": [], "thresholdI": [
    ], "lenKinIds": [], "sensitivityValue": []}
    resultsMatrix["fileName"].append(fileName)
    resultsMatrix["goodnessValue"].append(goodnessValue)
    resultsMatrix["thresholdI"].append(thresholdIndex)
    resultsMatrix["lenKinIds"].append(lenKinIds)
    resultsMatrix["sensitivityValue"].append(sensitiveValue)

    resultsMatrix.update(familyRangeIndexMatrix)

    return resultsMatrix

    # resultsMatrix.update(familyRangeIndexMatrix)
    # resultsMatrix = pd.DataFrame(resultsMatrix)

    # with open('file.txt', 'w') as file:
    #     # use `pickle.loads` to do the reverse
    #     file.write(str(resultsMatrix))
    #     file.close()
    # print(resultsMatrix)


if __name__ == "__main__":

    # __________MODIFY HERE THE INPUT WORKING FOLDER__________
    inputFolder = '../results/R06_hmmerTestCompareDataset/2019-05-28_001/'

    groupToAnalyze = 'STE'

    # analyzeHmmOut = ()

    # output_hmm = '../results/R06_hmmerTestCompareDataset/2019-05-28_001/msa_row4_col4.tblout'
    sequences_recognized = "our_kin.fasta"

    resultsMatrix = {}
    final_list = []
    i = 0
    for root, dirs, files in os.walk(inputFolder):
        for file in files:
            if (file.endswith(".tblout")):
                pathInputFile = inputFolder + file
                resultsMatrix = (analyse_hmm_out(
                    pathInputFile, sequences_recognized, resultsMatrix))
                print(pathInputFile)
                resultsMatrix2 = pd.DataFrame(resultsMatrix).values
                l = resultsMatrix2.shape[1]
                final_list.append(resultsMatrix2.reshape((l,)))
                print(i)
                i += 1
    print(final_list)

    with open('../results/R07_ClassifierResults/resultsMatrix.txt', 'w') as file:
        #     # use `pickle.loads` to do the reverse
        for element in resultsMatrix.keys():
            file.write(str(element) + " ")
        file.write("\n")
        for line in final_list:
            for element in line:
                file.write(str(element) + " ")
            file.write("\n")
            # file.write(str(final_list))
        file.close()

    # resultsMatrix = pd.DataFrame(resultsMatrix)
    # print(resultsMatrix)
# print(df)


# %%
