from Bio import SeqIO

kin = list(SeqIO.parse("kinase_dataset.fasta", "fasta"))

kin_descriptions = [kin.description.split()[1:] for kin in kin]

#kin_dict = {}

#for k in kin_descriptions:
#    kin_dict[k[0]] = (k[1],k[2])
    
out = open("kin_id_lab.txt","w")
out.write("Id Family Subfamily\n")

for k in kin_descriptions:
    out.write(k[0] + " " + k[1] + " " + k[2] + "\n")
out.close()


