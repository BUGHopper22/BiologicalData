from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np

alignment = AlignIO.read(open("msa200.clw"), "clustal")
# l è l'array costruito con ogni riga che corrisponde ad una sequenza
l = np.array([np.array(list(str(s.seq))) for s in alignment])
# l=="-" crea una matrice booleana con 1 ogni volta che trova un "-", sum somma tutti gli uni per ogni colonna, e facendo /l.shape[0] trovo la percentuale di - sul totale delle righe. 
# c sarà una matrice ad un vettore riga contenete le percentuali di "-" trovati per ogni colonna
c = np.sum(l == "-", axis = 0)/l.shape[0]

threshold = 0.1

# creo un array contenente i numeri da zero fino al numero di colonne-1
idx = np.arange(c.shape[0])
# 
idx = idx[c < threshold]

new_l = l[:,idx]

seq_n = [None]*len(alignment)

for i in range(len(alignment)):
    s = list(new_l[i,:])
    s = [str(el) for el in s]
    seq_n[i] = "".join(s)

msa = [None]*len(alignment)

for i in range(len(alignment)):
    msa[i] = SeqRecord(seq = Seq(seq_n[i],SingleLetterAlphabet()), id = alignment[i].id, 
       name = alignment[i].name, description = alignment[i].description, dbxrefs = alignment[i].dbxrefs)

align = MultipleSeqAlignment(msa)

out = open("align_cleaned.clw", "w")
AlignIO.write(align, out, "clustal")