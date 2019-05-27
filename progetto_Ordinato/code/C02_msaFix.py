
# input: R03_msa200.clw contains the multiple sequence allignement of all sequences
# output: R04_msa200Cleaned.clw
# function: remove row and column with too many gaps

from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
#import matplotlib.pyplot as plt


def cut_alignment(alignment, threshold_c=0, threshold_r=0):
    """
        the variable threshold_c is the maximum number of dashes accepted in the columns
        the variable threshold_r is the maximum number of dashes accepted in the rows
    """
    l = np.array([np.array(list(str(s.seq))) for s in alignment])

    c = np.sum(l == "-", axis=0)/l.shape[0]

    if threshold_c == 0:
        threshold_c = 0.1

    idx = np.arange(c.shape[0])
    idx = idx[c < threshold_c]

    #c_x = c[idx]
    # plt.plot(idx,c_x)
    new_l = l[:, idx]

    r = np.sum(new_l == "-", axis=1)/new_l.shape[1]

    if threshold_r == 0:
        threshold_r = 0.06

    idx = np.arange(r.shape[0])
    idx = idx[r > threshold_r]

    new_l[idx, :] = '?'  # check to not include undesired sequences

    seq_n = [None]*len(alignment)

    for i in range(len(alignment)):
        if i not in idx:
            s = list(new_l[i, :])
            s = [str(el) for el in s]
            seq_n[i] = "".join(s)

    msa_tmp = [None]*len(alignment)
    msa = []

    for i in range(len(alignment)):
        if i not in idx:
            msa_tmp[i] = SeqRecord(seq=Seq(seq_n[i], SingleLetterAlphabet()), id=alignment[i].id,
                                   name=alignment[i].name, description=alignment[i].description, dbxrefs=alignment[i].dbxrefs)
            msa.append(msa_tmp[i])

    
    align = MultipleSeqAlignment(msa)
    return align


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-ai", "--alignment_input", required=True,
                    help="alignment input file")
    ap.add_argument("-ao", "--alignment_output", required=True,
                    help="alignment output file")
    ap.add_argument("-tc", "--threshold_column", required=False,
                    help="maximum number of dashes accepted in the columns")
    ap.add_argument("-tr", "--threshold_row", required=True,
                    help="maximum number of dashes accepted in the rows")
    args = vars(ap.parse_args())

    alignment_input_file = args["alignment_input"]
    alignment_output_file = args["alignment_output"]
    threshold_c = float(args["threshold_column"])
    threshold_r = float(args["threshold_row"])

    alignment_input = AlignIO.read(open(alignment_input_file), "clustal")

    alignment_output = cut_alignment(alignment_input, threshold_c, threshold_r)

    outfile = open(alignment_output_file, "w")
    treshold = open('threshold.txt',"w")
    result.write("threshold_column: {}\n".format(threshold_c))
    result.write("threshold_row: {}\n".format(threshold_r))


    
    AlignIO.write(alignment_output, outfile, "clustal")
    outfile.close()
