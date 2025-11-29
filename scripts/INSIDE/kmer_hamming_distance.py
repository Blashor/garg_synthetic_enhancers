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


def kmer_hamming(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def kmer_count(seq, seq0, kmer_to_number, window=2):
    kmer_count_arr = np.zeros(len(kmer_to_number))
    kmer_hamming_arr = [[] for _ in range(len(kmer_to_number))]
    for i in range(len(seq) - window):
        kmer = seq[i : window + i]
        kmer0 = seq0[i : window + i]
        ham_dist = kmer_hamming(kmer0, kmer)
        kmer_r = rc(kmer)
        if kmer_to_number[kmer] < kmer_to_number[kmer_r]:
            kmer_count_arr[kmer_to_number[kmer]] += 1
            kmer_hamming_arr[kmer_to_number[kmer]].append(ham_dist)
        else:
            kmer_count_arr[kmer_to_number[kmer_r]] += 1
            kmer_hamming_arr[kmer_to_number[kmer_r]].append(ham_dist)
    return kmer_count_arr, kmer_hamming_arr


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
window = 8  # change this to desired size

kmer_labels = ["".join(p) for p in product(letters, repeat=window)]
kmer_to_number = {kmer: i for i, kmer in enumerate(kmer_labels)}


gen_sums = {}
kmer_hamm_total = [[] for _ in range(len(kmer_to_number))]
gen = "G90"
for line_rep in line_seq_holder:
    seq0 = line_seq_holder[line_rep]["G1"][0][0]
    seq = line_seq_holder[line_rep][gen][0][0]
    kmer_list, kmer_hamming_arr = kmer_count(seq, seq0, kmer_to_number, window=window)
    kmer_hamm_total = [a + b for a, b in zip(kmer_hamm_total, kmer_hamming_arr)]
    if gen in gen_sums:
        gen_sums[gen] += kmer_list
    else:
        gen_sums[gen] = np.array(kmer_list)


gens_to_plot = ["G90"]

# Convert labels to numpy array if not already
labels = np.array(kmer_labels)

violins = []
v_labels = []
for gen in gens_to_plot:
    y_counts = gen_sums[gen]
    threshold = np.percentile(y_counts, 99.9)

    for i, label in enumerate(labels):
        if y_counts[i] >= threshold:
            violins.append(kmer_hamm_total[i])
            v_labels.append(label)
    import pandas as pd

    # Flatten data into long form
    data = []
    for label, values in zip(v_labels, violins):
        for val in values:
            data.append({"kmer": label, "hamming": val})

    df = pd.DataFrame(data)

    plt.figure(figsize=(12, 6))
    sns.boxplot(x="kmer", y="hamming", data=df)
    plt.xticks(rotation=90)
    plt.title("Top 1% k-mers Hamming distribution")
    plt.show()
