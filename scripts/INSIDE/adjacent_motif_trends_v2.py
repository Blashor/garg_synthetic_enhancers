gen_add = 0
surv = 4
line_seq_holder = {}
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math


def find_differences(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Strings must have the same length")

    differences = []
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            differences.append(i)

    return differences


def find_consensus_diffs(diff_holder, str1):
    differences = []
    for i, base in enumerate(str1):
        if base not in diff_holder[i]:
            diff_holder[i].append(base)
            differences.append(i)
    return differences


with open("INSIDE_jun13.txt") as file:
    current_stat = 0
    gen_tracker = {}
    replicates_mean = {}
    itr = 0
    first_lines = True
    for line in file:
        line = line.strip()
        if "NEW_FILE" in line:
            gen_add += 1

        if "" not in line:
            continue

        if "Gen:" not in line:
            continue

        cell = line.split("\t")
        cell[0] = "Gen: " + str(int(cell[0].split("Gen:")[-1]) + gen_add)
        current_stat += float(cell[3])
        if cell[1] not in line_seq_holder:
            line_seq_holder[cell[1]] = {}
        if cell[0] not in line_seq_holder[cell[1]]:
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
snp_dist = []
score_dist = []
for line in line_seq_holder:
    last_seqs = []
    heatmap = []
    # print("NEW LINE")
    last_gen_score = 0
    last_gen_scores = [0, 0, 0, 0]
    score_strength = 0
    diff_holder = []
    for seq_pos in range(1024):
        diff_holder.append([])
    last_gen_snps = []

    for gen in line_seq_holder[line]:
        current_seqs = line_seq_holder[line][gen]
        # score_strength = math.log2(float(current_seqs[0][1]) / last_gen_score)

        # score_strength = float(current_seqs[0][1]) - last_gen_score
        print(gen, score_strength)
        current_diff_across = [0] * len(current_seqs[0][0])
        compound_snps = []
        last_gen_adders = []
        for i, c_seq in enumerate(current_seqs):
            snp_poses = find_consensus_diffs(diff_holder, c_seq[0])
            print(snp_poses)
            print(len(snp_poses))
            lgs = float(last_gen_scores[i])
            c_score = math.sqrt(float(current_seqs[i][1]))
            last_gen_adders.append(c_score)
            score_strength = c_score - lgs
            for pos in snp_poses:
                for last_pos in last_gen_snps:
                    if last_pos < 512 and pos < 512:
                        snp_dist.append((last_pos - pos))
                        score_dist.append(score_strength)
                    if last_pos > 512 and pos > 512:
                        snp_dist.append((last_pos - pos))
                        score_dist.append((score_strength))
                compound_snps.append(pos)
                current_diff_across[pos] += 1
        last_gen_snps = compound_snps
        last_gen_scores = last_gen_adders

        heatmap.append(current_diff_across)

        last_seqs = current_seqs
plt.figure()

sns.histplot(data=heatmap)
plt.show()
import pandas as pd

# Assume you already have x and y defined
x = snp_dist
y = score_dist

# Bin x into 64 bp intervals from -512 to 512
bin_edges = np.arange(0, 513, 4)
x_binned = pd.cut(x, bins=bin_edges)

# Create DataFrame
df = pd.DataFrame({"x_bin": x_binned, "y": y}).dropna()

# Plot
sns.boxplot(data=df, x="x_bin", y="y")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
# plt.scatter(snp_dist, score_dist, s=2, alpha=0.2)
sns.kdeplot(x=snp_dist, y=score_dist, s=2, alpha=0.2)
plt.show()

# Where did mutations occur
# delta odds ratio by position
