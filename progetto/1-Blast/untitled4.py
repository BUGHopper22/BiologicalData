import numpy as np
import pandas as pd

msa_file = "MSAUniref90.fasta"

strings = [line.split("\t") for line in open(msa_file) if len(line) > 0][1:-1]

keys = []
scores = []
alignments = []

for string in strings:
    s = string[0].split(" ")
    k = s[0]
    sc = string[-1]
    al = s[-1]
    if len(k) > 0 and k[0] == 'U':
        keys.append(k)
        scores.append(sc)
        alignments.append(al)

scores = np.array(scores)
alignments = np.array(alignments)
unique_keys = np.unique(keys)

msa = {}

for i in range(len(unique_keys)):
   al = alignments[np.arange(0,len(keys),len(unique_keys)) + i]
   s = ""
   for j in range(len(al)):
       s = s + al[j]
   msa[unique_keys[i]] = list(s)
   
msa_matrix = pd.DataFrame(msa).values.T

pd.DataFrame(msa_matrix,index = msa.keys()).to_csv("alignment.txt",sep = ' ')