from Bio import SearchIO

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
    
    