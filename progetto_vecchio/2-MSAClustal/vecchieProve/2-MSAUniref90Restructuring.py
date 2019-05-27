import numpy as np
import pandas as pd

msa_file = "MSAUniref90.fasta"
# crea un vettore colonna contenete stringhe (con id sequenza e score da dividere successivamente)
strings = [line.split("\t") for line in open(msa_file) if len(line) > 0][1:-1]

keys = []
scores = []
alignments = []

for string in strings:
    #split è una funzione che separa la stringa in n colonne, n è un numero deciso in base al carattere di separazion (in questo caso " ")
    # in questo caso ha creato un array[0,1,0] dove s[0]=id,s[1]=seq,s[2]=score?
    s = string[0].split(" ") 
    k = s[0]
    sc = string[-1]
    al = s[-1]
    if len(k) > 0 and k[0] == 'U':
        keys.append(k)
        scores.append(sc)
        alignments.append(al)
# questo sarebbe un cast da array normale a numpyarray?
scores = np.array(scores)
alignments = np.array(alignments)
unique_keys = np.unique(keys) #unique è un array?

#dizionario con chiave e valore
msa = {}

for i in range(len(unique_keys)):
   al = alignments[np.arange(0,len(keys),len(unique_keys)) + i]
   s = ""
   for j in range(len(al)):
       s = s + al[j]
   msa[unique_keys[i]] = list(s)
   
msa_matrix = pd.DataFrame(msa).values.T

pd.DataFrame(msa_matrix,index = msa.keys()).to_csv("alignment.txt",sep = ' ')