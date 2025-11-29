import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import math
import matplotlib.cm as cm
import random


# Example usage


survivors = 4
bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
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
# print(line_seq_holder)
with open("fimo_40/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, surv_rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            gen_n = int(gen.split("G")[-1])
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

violin_dict = {}
color_motif = {}
clrs = sns.color_palette("deep", 10) + sns.color_palette("husl", 191)
clr_i = 0
graphing_obj = {}
node_scores = {}
generations = {}
itr = 0
motif_appears_in_replicates = {}
motif_sets = []
allowed_motifs = set(
    [
        "BACH2",
        "NFE2",
        "SOX2",
        "FOXC1",
        "RFX6",
        "FOXD3",
        "FOXA2",
        "SOX10",
        "MAFB",
        "RFX3",
        "ZIC3",
        "PO2F2",
        "MAFK",
        "ZBT17",
        "PO5F1",
        "NANOG",
        "LEF1",
        "PO3F1",
        "SOX3",
        "MAFG",
        "NF2L2",
        "SOX9",
        "BACH1",
        "PO2F1",
        "RFX2",
        "ZIC2",
        "RFX1",
        "NR5A2",
        "SOX4",
        "MAFF",
        "CTCF",
        "CTCFL",
    ]
)

pop_level_data = {}
motif_set = {}
rl_count = {}
for line in line_seq_holder:
    old_motifs = set()
    for gen in line_seq_holder[line]:
        gen_num = int(gen.replace("G", ""))
        last_gen = f"G{gen_num-1}"
        if gen_num != 0:
            try:
                itr += 1
                # print(gen)
                score = math.sqrt(float(line_seq_holder[line][gen][0][1]))
                last_score = math.sqrt(float(line_seq_holder[line][last_gen][0][1]))
                seq = line_seq_holder[line][gen][0][0]
                motifs = line_seq_holder[line][gen][1:]
                # print(motifs)
                right = set()
                left = set()
                # print(gen, score, motifs)
                motif_count = {}
                for motif, m_pos in motifs:
                    if motif not in motif_count:
                        motif_count[motif] = 0
                    motif_count[motif] += 1
                    motif_set[motif] = []
                    if m_pos <= 512:
                        # left
                        left.add(motif)
                    elif m_pos > 512:
                        right.add(motif)
                if motif_count["SOX2"] == 2:
                    side = ""
                    if "SOX2" in right:
                        side += "r"
                    if "SOX2" in left:
                        side += "l"
                    if side not in rl_count:
                        rl_count[side] = []
                    rl_count[side].append(score)
                    print(side, score)

                for m in right:
                    if motif_count[m] != 2:
                        right.remove(m)
                for m in left:
                    if motif_count[m] != 2:
                        left.remove(m)
                both = right.intersection(left)
                right_only = right - left
                left_only = left - right
                if gen not in pop_level_data:
                    pop_level_data[gen] = {}
                for m in both:
                    if m not in pop_level_data[gen]:
                        pop_level_data[gen][m] = {"both": 0, "right": 0, "left": 0}
                    pop_level_data[gen][m]["both"] += 1
                for m in right_only:
                    if m not in pop_level_data[gen]:
                        pop_level_data[gen][m] = {"both": 0, "right": 0, "left": 0}
                    pop_level_data[gen][m]["right"] += 1
                for m in left_only:
                    if m not in pop_level_data[gen]:
                        pop_level_data[gen][m] = {"both": 0, "right": 0, "left": 0}
                    pop_level_data[gen][m]["left"] += 1

                # print(itr, motifs)

            except Exception:
                pass

for rl in rl_count:
    print(rl, len(rl_count[rl]), np.median(rl_count[rl]))
