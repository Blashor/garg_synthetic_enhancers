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
dir20 = "/Users/blake/Downloads/INSIDE20b"

in20s = os.listdir(dir20)
for f20 in in20s:
    with open(f"{dir20}/{f20}") as file:
        for line in file:
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

                cell[2] = cell[2].split("'")[1] + cell[2].split("'")[3]
                line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
letters = ["A", "C", "T", "G"]
window = 7  # change this to desired size

kmer_labels = ["".join(p) for p in product(letters, repeat=window)]
kmer_to_number = {kmer: i for i, kmer in enumerate(kmer_labels)}


gen_sums = {}

gen_ors = []
np_calcs = []
all_lines = []

for i, line_rep in enumerate(line_seq_holder):
    line_or = []
    seq_len = len(line_seq_holder[line_rep][f"G1"][0][0])
    np_calcs.append(seq_len)

    for num in range(50):
        # seq = line_seq_holder[line_rep][gen][0][0]
        e_or = math.sqrt(float(line_seq_holder[line_rep][f"G{num+1}"][0][1])) / math.sqrt(
            float(line_seq_holder[line_rep][f"G{1}"][0][1])
        )
        e_or = math.sqrt(float(line_seq_holder[line_rep][f"G{num+1}"][0][1]))
        line_or.append((e_or))
    all_lines.append(line_or)
    plt.plot(list(range(50)), line_or, color="black", alpha=0.03)
all_lines = np.array(all_lines)  # shape: (num_lines, 50)
percentiles = np.percentile(all_lines, [25, 50, 75], axis=0)  # shape: (3, 50)

# ---- Plot the percentile curves ----
gens = np.arange(50)
plt.plot(gens, percentiles[1], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[0], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[2], color="white", label="Median", linewidth=2, linestyle="--")
# plt.fill_between(gens, percentiles[0], percentiles[2], color="red", alpha=0.2, label="25–75% range")

plt.show()
print(np.percentile(np.array(np_calcs), [25, 50, 100]))

gen_ors = []
np_calcs = []
all_lines = []

for i, line_rep in enumerate(line_seq_holder):
    line_or = []
    seq_len = len(line_seq_holder[line_rep][f"G1"][0][0])
    np_calcs.append(seq_len)
    for num in range(50):
        # seq = line_seq_holder[line_rep][gen][0][0]
        e_or = math.sqrt(float(line_seq_holder[line_rep][f"G{num+1}"][0][1])) / math.sqrt(
            float(line_seq_holder[line_rep][f"G{1}"][0][1])
        )
        # e_or = math.sqrt(float(line_seq_holder[line_rep][f"G{num+1}"][0][1]))
        line_or.append((e_or))
    all_lines.append(line_or)
    plt.plot(list(range(50)), line_or, color="black", alpha=0.005)
all_lines = np.array(all_lines)  # shape: (num_lines, 50)
percentiles = np.percentile(all_lines, [25, 50, 75], axis=0)  # shape: (3, 50)

# ---- Plot the percentile curves ----
gens = np.arange(50)
plt.yscale("log")
plt.plot(gens, percentiles[1], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[0], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[2], color="white", label="Median", linewidth=2, linestyle="--")
plt.show()
print(np.percentile(np.array(np_calcs), [25, 50, 100]))

gen_ors = []
np_calcs = []
all_lines = []

for i, line_rep in enumerate(line_seq_holder):
    line_or = []
    seq_len = len(line_seq_holder[line_rep][f"G1"][0][0])
    np_calcs.append(seq_len)
    for num in range(50):
        seq = line_seq_holder[line_rep][f"G{num+1}"][0][0]
        # e_or = math.sqrt(float(line_seq_holder[line_rep][f"G{num+1}"][0][1]))
        count = (seq.count("G") + seq.count("C")) / len(seq)
        line_or.append(count)
    all_lines.append(line_or)
    plt.plot(list(range(50)), line_or, color="black", alpha=0.01)
all_lines = np.array(all_lines)  # shape: (num_lines, 50)
percentiles = np.percentile(all_lines, [25, 50, 75], axis=0)  # shape: (3, 50)

# ---- Plot the percentile curves ----
gens = np.arange(50)
plt.plot(gens, percentiles[1], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[0], color="white", label="Median", linewidth=2, linestyle="--")
plt.plot(gens, percentiles[2], color="white", label="Median", linewidth=2, linestyle="--")
plt.show()
print(np.percentile(np.array(np_calcs), [25, 50, 100]))

# xgen_ors.append(e_or)
