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


with open("INSIDE_gen80_rep40.out") as file:
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

bg_or_lines = []
bg_snp_lines = []

for sample_i, line in enumerate(line_seq_holder):
    heatmap = []
    or_holder = [0] * 1024
    snp_num_holder = [0] * 1024
    for itr, gen in enumerate(line_seq_holder[line]):
        print(gen)
        if itr != 0 and itr < 149:
            last_gen = f"Gen: {150}"
            snp_holder = []
            old_scores = []
            new_scores = []
            new_diffs = []
            current_diff_across = [0] * len(line_seq_holder[line][gen][0][0])
            for seq_pos in range(1024):
                snp_holder.append([])
            for seq, score in line_seq_holder[line][last_gen]:
                old_scores.append(float(score))
                difs = find_consensus_diffs(snp_holder, seq)
            for seq, score in line_seq_holder[line][gen]:
                new_scores.append(float(score))
                snp_difs = find_consensus_diffs(snp_holder, seq)
                # print(score, difs)
                score_delta = math.sqrt(float(score)) - math.sqrt(np.median(old_scores))
                for pos in snp_difs:
                    or_holder[pos] += score_delta
                    snp_num_holder[pos] += 1
                    new_diffs.append(pos)
                    # current_diff_across[pos] += -1

                for pos, snp in enumerate(seq):
                    if snp == "A":
                        current_diff_across[pos] = 1
                    elif snp == "G":
                        current_diff_across[pos] = 2
                    elif snp == "C":
                        current_diff_across[pos] = 3
                    elif snp == "T":
                        current_diff_across[pos] = 4
            # score_delta = math.sqrt(np.median(new_scores)) - math.sqrt(np.median(old_scores))

            heatmap.append(current_diff_across)
    # print(or_holder)
    # print(snp_num_holder)
    chrom = f"chr{sample_i+1}"
    for pos_i, score in enumerate(or_holder):
        bg_or_lines.append(f"{chrom}\t{pos_i}\t{pos_i+1}\t{score}\n")
    for pos_i, score in enumerate(snp_num_holder):
        bg_snp_lines.append(f"{chrom}\t{pos_i}\t{pos_i+1}\t{score}\n")

    plt.figure()
    sns.heatmap(heatmap, cmap="magma")

    #
    plt.figure()
    plt.plot(range(len(or_holder)), or_holder)
    plt.plot(range(len(or_holder)), snp_num_holder, linewidth=0.5, color="red", alpha=0.2)
    plt.show()

with open("igvish_snp_num.bedgraph", "w") as file:
    file.writelines(bg_snp_lines)
with open("igvish_or_delta.bedgraph", "w") as file:
    file.writelines(bg_or_lines)
# Where did mutations occur
# delta odds ratio by position
