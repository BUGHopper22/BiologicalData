from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

alignment = AlignIO.read(open("msa200.clw"), "clustal")

l = np.array([np.array(list(str(s.seq))) for s in alignment])

c = np.sum(l == "-", axis = 0)/l.shape[0]

threshold_c = 0.1

idx = np.arange(c.shape[0])
idx = idx[c < threshold_c]

c_x = c[idx]

plt.plot(idx,c_x)
new_l = l[:,idx]

r = np.sum(new_l == "-", axis = 1)/new_l.shape[1]

threshold_r = 0.06

idx = np.arange(r.shape[0])
idx = idx[r > threshold_r]

new_l[idx,:] = '?' # check to not include undesired sequences

seq_n = [None]*len(alignment)

for i in range(len(alignment)):
    if i not in idx:
        s = list(new_l[i,:])
        s = [str(el) for el in s]
        seq_n[i] = "".join(s)

msa_tmp = [None]*len(alignment)
msa = []

for i in range(len(alignment)):
    if i not in idx:
        msa_tmp[i] = SeqRecord(seq = Seq(seq_n[i],SingleLetterAlphabet()), id = alignment[i].id, 
         name = alignment[i].name, description = alignment[i].description, dbxrefs = alignment[i].dbxrefs)
        msa.append(msa_tmp[i])

    
align = MultipleSeqAlignment(msa)

out = open("align_cleaned.clw", "w")
AlignIO.write(align, out, "clustal")