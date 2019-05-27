from Bio import SeqIO

cl_file = "kinase_cd_hit.txt"

d = [line for line in open(cl_file)]

i = 0
k = 0
cluster_ids = {}

while i < len(d):
    if d[i][0] == '>':
        k = k + 1
        cluster_ids[k] = []
    else:
        cluster_ids[k].append(d[i][d[i].find(">") + 1 : d[i].find("...")])
    i = i + 1
    
ste_kin = list(SeqIO.parse("ste_kin.fasta", "fasta"))

ste_kin_names = [kin.name for kin in ste_kin]

for key in cluster_ids:
    cluster = []
    for name in cluster_ids[key]:
        cluster.append(ste_kin[ste_kin_names.index(name)])
    SeqIO.write(cluster, str(key) + "_cluster.fasta", "fasta")
        