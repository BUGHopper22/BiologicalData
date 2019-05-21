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