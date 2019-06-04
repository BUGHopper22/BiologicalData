# %%
from Bio import SearchIO
import numpy as np
from sklearn.metrics import confusion_matrix
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#from Bio import SeqIO
import pandas as pd
import os
#import _pickle as pickle


def insertDataInArrays(qresults, scores, group, family):
    for i, qresult in enumerate(qresults.hits):
        scores[i] = qresult.bitscore
        desc = qresult.description.split()
        group[i] = desc[1]
        family[i] = desc[2]

# create three arrays with statistical value needed


def buildPerformanceMeasuresLists(scores, myFamilyTrue):

    maxBitScore = int(np.max(scores)+10)
    steps = np.min(-np.diff(scores)) - 0.0001
    if steps <= 0:
        steps = 0.09
        
    length = len(np.arange(0,maxBitScore,steps))
    
    accuracies = [0]*length
    sensitivities = [0]*length
    specificities = [0]*length
    thresholds = [0]*length

    for i, j in enumerate(np.arange(0, maxBitScore,steps)):

        myFamilyPred = np.int32(scores > j)

        cm = confusion_matrix(myFamilyTrue, myFamilyPred)

        tn = cm[0, 0]
        fp = cm[0, 1]
        fn = cm[1, 0]
        tp = cm[1, 1]

        acc = (tp + tn) / (tp + fp + fn + tn)
        sens = tp / (tp + fn)
        spec = tn / (tn + fp)
        
        thresholds[i] = j
        accuracies[i] = acc
        sensitivities[i] = sens
        specificities[i] = spec

    return accuracies, sensitivities, specificities, thresholds


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
    
    myGroupUniqueFamilies = np.unique(family[myGroupBoolArray])

    # true become 1, false become 0
    myGroupTrue = np.int32(myGroupBoolArray)

    # create three list with statistical value needed
    accuracies, sensitivities, specificities, thresholds = buildPerformanceMeasuresLists(
        scores, myGroupTrue)
    # print(accuracies)

    thresholdIndex = np.argmax(accuracies)
    thresholdValue = thresholds[thresholdIndex]
    print(thresholdIndex)
    # print('thresholdIndex')
    # print(thresholdIndex)

    kin_ids = []
    kin_ids = insertDataInKin_ids(kin_ids, thresholdValue, qresults)

    fileName = out_hmm
    goodnessValue = np.sum((myGroupIndexies + 1)**2)

    kinRec = len(kin_ids)
    sensitivityValue = sensitivities[thresholdIndex]
    # print('sensitiveValue')
    # print(sensitiveValue)
    
    familiesCaptured = family[:kinRec]
    
    resultsMatrix = {}
    resultsMatrix["fileName"] = [fileName]
    resultsMatrix["goodnessValue"] = [goodnessValue]
    resultsMatrix["sensitivityValue"] = [sensitivityValue]
    resultsMatrix["threshold"] = [thresholdValue]
    resultsMatrix["kinRec"] = [kinRec]
    for fam in myGroupUniqueFamilies:
        indices_fam = np.arange(N)[family == fam]
        resultsMatrix[fam + "-sens"] = [np.sum(familiesCaptured == fam)/len(indices_fam)]
        resultsMatrix[fam + "-BI"] = [indices_fam[0]]
        resultsMatrix[fam + "-WI"] = [indices_fam[-1]]
 
    print(resultsMatrix)
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
        print(files)
        for file in files:
            if (file.endswith(".tblout")):
                pathInputFile = inputFolder + file
                resultsMatrix = (analyse_hmm_out(
                    pathInputFile, sequences_recognized, resultsMatrix))
                print(pathInputFile)
                resultsMatrix2 = pd.DataFrame(resultsMatrix).values
                print(resultsMatrix2)
                l = resultsMatrix2.shape[1]
                final_list.append(resultsMatrix2.reshape((l,)))
                # print('resultsMatrix')
                # print(resultsMatrix)
                # print(type(resultsMatrix))

                # resultsMatrix = pd.DataFrame(resultsMatrix)
                # print('resultsMatrix')
                # print(resultsMatrix)
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
