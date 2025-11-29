bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from itertools import product


def rc(dna_seq):
    # Create complement mapping
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    # Get complement and reverse
    return dna_seq.translate(complement)[::-1]


def kmer_count(seq, kmer_to_number, window=2):
    kmer_count_arr = np.zeros(len(kmer_to_number))
    for i in range(len(seq) - window):
        kmer = seq[i : window + i]
        kmer_r = rc(kmer)
        if kmer_to_number[kmer] < kmer_to_number[kmer_r]:
            kmer_count_arr[kmer_to_number[kmer]] += 1
        else:
            kmer_count_arr[kmer_to_number[kmer_r]] += 1
    return kmer_count_arr


tsv_lines = []
# INSIDE_256_g117.out
with open("INSIDE_jun13.txt") as file:
    current_stat = 0
    gen_tracker = {}
    replicates_mean = {}
    itr = -1
    first_lines = True
    for line in file:
        line = line.strip()

        if "NEW_FILE" in line:
            gen_add += 1

        if "" not in line:
            continue

        if "Gen:" not in line:
            continue
        itr += 1
        surv_num = (itr % survivors) + 1
        cell = line.split("\t")
        cell[0] = "G" + str(int(cell[0].split("Gen:")[-1]) + gen_add)
        current_stat += float(cell[3])
        cell[1] += "_" + str(surv_num)
        if cell[1] not in line_seq_holder:
            line_seq_holder[cell[1]] = {}
        if cell[0] not in line_seq_holder[cell[1]]:
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
letters = ["A", "C", "T", "G"]
window = 5  # change this to desired size

kmer_labels = ["".join(p) for p in product(letters, repeat=window)]
kmer_to_number = {kmer: i for i, kmer in enumerate(kmer_labels)}


gen_sums = {}

for line_rep in line_seq_holder:
    for gen in line_seq_holder[line_rep]:
        seq = line_seq_holder[line_rep][gen][0][0]
        kmer_list = kmer_count(seq, kmer_to_number, window=window)

        if gen in gen_sums:
            gen_sums[gen] += kmer_list
        else:
            gen_sums[gen] = np.array(kmer_list)
top_n = 30  # number of top k-mers to display

gens_to_plot = ["G30", "G60", "G90", "G120"]

# Convert labels to numpy array if not already
labels = np.array(kmer_labels)
g1_counts = gen_sums["G30"]

for gen in gens_to_plot:
    y_counts = gen_sums[gen]

    plt.figure(figsize=(6, 4))
    plt.scatter(g1_counts, y_counts)
    # plt.plot(g1_counts, y_counts, linestyle="--", color="gray")  # optional line
    plt.title(f"{gen} vs G1 K-mers")
    plt.xlabel("G1 Counts")
    plt.ylabel(f"{gen} Counts")

    # Annotate points with k-mer labels
    # Determine top 1% of y_counts
    threshold = np.percentile(y_counts, 98)

    # Annotate only top 1%
    for i, label in enumerate(labels):
        if y_counts[i] >= threshold:
            plt.text(g1_counts[i], y_counts[i], label, fontsize=8, ha="right", va="bottom")

    plt.show()
