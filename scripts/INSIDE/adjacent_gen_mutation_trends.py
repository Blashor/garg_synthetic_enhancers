gen_add = 0
surv = 4
line_seq_holder = {}
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


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
control_dist = []
import random
import math

for line in line_seq_holder:
    last_seqs = []
    heatmap = []
    # print("NEW LINE")
    last_gen_score = 1
    score_strength = 0
    diff_holder = []
    for seq_pos in range(1024):
        diff_holder.append([])
    last_gen_snps = []

    for gen in line_seq_holder[line]:
        current_seqs = line_seq_holder[line][gen]
        score_strength = math.sqrt(float(current_seqs[0][1])) - last_gen_score
        print(gen, score_strength)
        last_gen_score = math.sqrt(float(current_seqs[0][1]))
        current_diff_across = [0] * len(current_seqs[0][0])
        compound_snps = []
        for c_seq in current_seqs:
            snp_poses = find_consensus_diffs(diff_holder, c_seq[0])
            print(snp_poses)
            print(len(snp_poses))
            if score_strength > 1.87:
                for pos in snp_poses:
                    for last_pos in last_gen_snps:
                        if last_pos < 512 and pos < 512:
                            snp_dist.append(last_pos - pos)
                            control_dist.append(last_pos - random.randint(0, 511))
                            score_dist.append(score_strength)
                        if last_pos > 512 and pos > 512:
                            snp_dist.append(last_pos - pos)
                            score_dist.append(score_strength)
                            control_dist.append(last_pos - random.randint(513, 1023))
                    compound_snps.append(pos)
                    current_diff_across[pos] += 1
        last_gen_snps = compound_snps

        heatmap.append(current_diff_across)

        last_seqs = current_seqs
print(len(snp_dist))


random_list = np.array([random.randint(0, 512) for _ in range(len(snp_dist))]) - np.array(
    [random.randint(0, 512) for _ in range(len(snp_dist))]
)
import numpy as np

print(np.percentile(score_dist, 99))
mean = np.mean(snp_dist)
std_dev = np.std(snp_dist)
lower, upper = -512, 512
target_size = len(snp_dist)

new_dist = []

samples = np.random.normal(loc=mean, scale=std_dev, size=target_size)
plt.figure()

sns.ecdfplot(data=control_dist, label="Ctrl", alpha=0.8)
sns.ecdfplot(data=samples, label="Normal", alpha=0.5)
sns.ecdfplot(data=snp_dist, label="SNP", alpha=0.8)
plt.legend()
plt.show()

# Where did mutations occur
# delta odds ratio by position
