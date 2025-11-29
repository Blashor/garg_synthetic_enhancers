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
np_heat = np.zeros([300, 1024])
for line in line_seq_holder:
    last_seqs = []
    heatmap = []
    # print("NEW LINE")
    last_gen_score = 1
    score_strength = 0
    diff_holder = []
    for seq_pos in range(1024):
        diff_holder.append([])
    current_diff_across = np.zeros([300, 1024])
    for gen in line_seq_holder[line]:
        gen_int = int(gen.split(" ")[-1])
        if gen_int < 300:
            current_seqs = line_seq_holder[line][gen]
            score_strength = float(current_seqs[0][1]) / last_gen_score
            print(gen, score_strength)
            last_gen_score = float(current_seqs[0][1])

            for c_seq in current_seqs:
                snp_poses = find_consensus_diffs(diff_holder, c_seq[0])
                print(len(snp_poses))
                for pos in snp_poses:
                    if pos > 511:
                        pos -= 512
                    current_diff_across[gen_int, pos] += 1

            # heatmap.append(current_diff_across)

            last_seqs = current_seqs
    np_heat += current_diff_across
plt.figure()
sns.heatmap(np_heat[2:, :512], cmap="magma")
x_tick_positions = np.arange(128, 512, 128)  # Adjust interval as needed
plt.xticks(ticks=x_tick_positions, labels=x_tick_positions, rotation=0, ha="center")
y_tick_positions = np.arange(30, 300, 30)  # Adjust interval as needed
plt.yticks(ticks=y_tick_positions, labels=y_tick_positions, rotation=0)
plt.xlabel("Position (bp)")
plt.ylabel("Generation")
plt.savefig("map_snps.png", dpi=300)
plt.show()

# Where did mutations occur
# delta odds ratio by position
