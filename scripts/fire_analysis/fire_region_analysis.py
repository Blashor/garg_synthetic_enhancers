from pyfaidx import Fasta
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import mpl_style
import numpy as np


#
# bed element gc content
#
def get_gc(seq):
    gc = 0
    if "N" not in seq:
        for s in seq:
            if s == "G" or s == "C":
                gc += 1
        return gc / len(seq)
    else:
        return None


def calculate_gc_dist(file_name):
    fasta = Fasta("/Users/blake/Downloads/GRCm38.primary_assembly.genome.fa")
    gc_dist = []

    with open(file_name) as file:
        for line in file:
            cell = line.split("\t")
            # print(cell[0], cell[1], cell[2])
            if "#chrom" not in cell[0]:
                try:
                    seq = str(fasta[cell[0]][int(cell[1]) : int(cell[2])])
                    gc_per = get_gc(seq)
                    if gc_per != None:
                        gc_dist.append(gc_per)
                except Exception as e:
                    pass

    return gc_dist


"""
os.system(
    f"bedtools shuffle -i ../data/3runs-FDR-FIRE-peaks.bed -g /Users/blake/Downloads/GRCm38.primary_assembly.genome.fa.fai > ../data/Fire_shuffle.bed"
)
"""
mpl_style.set_style(figsize=(4.5, 4.5))
gc_dist = calculate_gc_dist("../data/3runs-FDR-FIRE-peaks.bed")
sns.ecdfplot(data=gc_dist, color=(157 / 255, 30 / 255, 30 / 255), label=f"FIRE, Median: {round(np.median(gc_dist),3)}")

gc_dist = calculate_gc_dist("/Users/blake/Downloads/SRX10040677_dnase.05.bed")
sns.ecdfplot(data=gc_dist, label=f"DNAseI, Median: {round(np.median(gc_dist),3)}")
gc_dist = calculate_gc_dist("/Users/blake/Downloads/SRX23682289_atac2.05.bed")
sns.ecdfplot(data=gc_dist, label=f"ATAC-seq, Median: {round(np.median(gc_dist),3)}")
# gc_dist_control = calculate_gc_dist("../data/Fire_shuffle.bed")
# sns.ecdfplot(data=gc_dist_control, color="black", label=f"Control, Median: {round(np.median(gc_dist_control),3)}")
plt.ylim([0, 1.1])
plt.xlim([0, 1.0])
plt.legend()
plt.savefig("../figures/2b_seq_gc_analysis.svg")
plt.show()
