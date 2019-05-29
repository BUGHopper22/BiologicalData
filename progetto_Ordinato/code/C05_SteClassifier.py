# input: HMMER output
# output: kinase_hmmer_cleaned.tblout
# function: 

from Bio import SearchIO
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import SeqIO

out_hmm = 'kinase_hmmer_cleaned.tblout'

qresults = next(SearchIO.parse(out_hmm, 'hmmer3-tab'))

N = len(qresults.hits)
scores = [0]*N
group = [None]*N
family = [None]*N

for i,qresult in enumerate(qresults.hits):
    scores[i] = qresult.bitscore
    desc = qresult.description.split()
    group[i] = desc[1]
    family[i] = desc[2]

group = np.array(group)
scores = np.array(scores)
y_true = np.int32(group == 'STE')

accuracies = [0]*800
sensitivities = [0]*800
specificities = [0]*800

for i,j in enumerate(np.arange(0,800)):
    
    y_pred = np.int32(scores > j)
    
    cm = confusion_matrix(y_true,y_pred)
    
    tn = cm[0,0]
    fp = cm[0,1]
    fn = cm[1,0]
    tp = cm[1,1]
    acc = (tp+tn)/(tp+fp+fn+tn)
    sens = tp/(tp+fn)
    spec = tn/(tn+fp)
    accuracies[i] = acc
    sensitivities[i] = sens
    specificities[i] = spec

fig_n = 1
fig = plt.figure(fig_n)
fig_n = fig_n + 1
gs = gridspec.GridSpec(1,1)
ax0 = fig.add_subplot(gs[0,0])
ax0.plot(accuracies,'g')
ax0.plot(sensitivities,'b')
ax0.plot(specificities,'r')
ax0.legend(['Accuracies', 'Sensitivities', 'Specificities'])

y = np.array([0,1])
m = np.argmax(accuracies)
x = np.array([m,m])

ax0.plot(x,y,'k')

i = 0
sensitivities = np.array(sensitivities)
threshold = m
#threshold = len(sensitivities[sensitivities > 0.94])
kin_ids = []

while qresults.hits[i].bitscore > threshold:
    desc = qresults.hits[i].description
    kin_id = desc.split()[0]
    kin_ids.append(kin_id)
    i = i + 1

kinase_dataset_path = "kinase_dataset.fasta"

kinase_dataset = list(SeqIO.parse(open(kinase_dataset_path), "fasta"))

kin_our = [None]*len(kin_ids)
i = 0

for kinase in kinase_dataset:
    if kinase.description.split()[1] in kin_ids:
        kin_our[i] = kinase
        i = i + 1
    
with open("our_kin.fasta", "w") as output_handle:
    SeqIO.write(kin_our, output_handle, "fasta")