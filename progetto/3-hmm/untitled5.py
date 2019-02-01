import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

out_hmm = 'kinase_hmmer_cleaned.ali'

l = [line.split(" ") for line in open(out_hmm)][15:375]

tmp = [f for f in l if len(f)>0]

families = np.array([f[-2] for f in l])

y_true = families == 'STE'
accuracies = [0]*len(families)
sensitivities = [0]*len(families)
specificities = [0]*len(families)



for i in range(800):
    
    y_pred = np.int32(scores > i)
    


for i,index in enumerate(range(len(families))):

    y_pred = np.concatenate([np.ones((index,)),np.zeros((len(l)-index,))])
    
    cm = confusion_matrix(y_true,y_pred)
    
    tp = cm[0,0]
    fp = cm[0,1]
    fn = cm[1,0]
    tn = cm[1,1]
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
ax0.plot(scores,accuracies,'go',markersize = 2)
ax0.plot(scores,sensitivities,'bo',markersize = 2)
ax0.plot(scores,specificities,'ro',markersize = 2)
ax0.legend(['Accuracies', 'Sensitivities', 'Specificities'])