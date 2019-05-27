from Bio import SeqIO
import pandas as pd

kin = list(SeqIO.parse("sequences.fasta", "fasta"))

kin_ids = [k.id for k in kin]

kin_id_lab = pd.read_csv("kin_id_lab.txt", sep = " ")

for i, k_id in enumerate(kin_ids):
    if k_id[:3] != "UPI":
        kin_ids[i] = kin_ids[i][3: 3 + kin_ids[i][3:].find("|")]

kin_ids = pd.DataFrame({"Id": kin_ids})
kin_ids = kin_ids.merge(kin_id_lab,how = "left")