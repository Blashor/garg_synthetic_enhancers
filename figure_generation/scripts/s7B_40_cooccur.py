import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import math
import matplotlib.cm as cm
import random


def plot_motif_cooccurrence(motif_sets):
    # Get all unique motifs
    unique_motifs = set()
    for motif_set in motif_sets:
        unique_motifs.update(motif_set)

    # Create a list and a mapping from motif to index
    unique_motifs = list(unique_motifs)
    motif_to_index = {motif: idx for idx, motif in enumerate(unique_motifs)}

    # Number of unique motifs
    num_motifs = len(unique_motifs)

    # Initialize the co-occurrence matrix
    cooccurrence_matrix = np.zeros((num_motifs, num_motifs), dtype=int)

    # Update the co-occurrence matrix
    for motif_set in motif_sets:
        for motif1 in motif_set:
            for motif2 in motif_set:
                if motif1 != motif2:
                    idx1 = motif_to_index[motif1]
                    idx2 = motif_to_index[motif2]
                    cooccurrence_matrix[idx1, idx2] += 1

    # Plot the co-occurrence matrix
    plt.figure(figsize=(10, 8))
    # sns.heatmap(cooccurrence_matrix, , yticklabels=unique_motifs, cmap="Blues")
    cluster = sns.clustermap(
        cooccurrence_matrix,
        xticklabels=unique_motifs,
        yticklabels=unique_motifs,
        method="ward",
        cmap="Blues",
        annot=True,
        fmt=".3g",
    )
    plt.setp(cluster.ax_heatmap.yaxis.get_majorticklabels(), fontsize=8)
    plt.setp(cluster.ax_heatmap.xaxis.get_majorticklabels(), fontsize=8)
    # plt.title("Motif Co-occurrence Matrix")
    plt.savefig("../figures/s8c_cooccur.svg")
    plt.show()


# Example usage


survivors = 4
bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
tsv_lines = []
# INSIDE_256_g117.out
with open("../data/INSIDE_jun13.txt") as file:
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
with open("../data/INSIDE_40_definitive/fimo.tsv") as file:
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


for line in line_seq_holder:
    old_motifs = set()
    for gen in line_seq_holder[line]:
        gen_num = int(gen.replace("G", ""))
        last_gen = f"G{gen_num-1}"
        if gen_num == 300:
            try:
                itr += 1
                print(gen)
                score = math.sqrt(float(line_seq_holder[line][gen][0][1]))
                last_score = math.sqrt(float(line_seq_holder[line][last_gen][0][1]))
                seq = line_seq_holder[line][gen][0][0]
                motifs = line_seq_holder[line][gen][1:]
                # print(gen, score, motifs)
                motif_holder = {}
                for motif, m_pos in motifs:
                    if motif not in motif_holder:
                        motif_holder[motif] = 0
                    motif_holder[motif] += 1
                motifs = set()
                for motif in motif_holder:
                    if motif in allowed_motifs:
                        motif_id = f"{motif}"  # _{motif_holder[motif]}"

                        motifs.add(motif_id)
                motif_sets.append(motifs)
                print(itr, motifs)

            except Exception:
                pass
plot_motif_cooccurrence(motif_sets)


# net_x(graphing_obj, node_scores, motif_sizes, generations, motif_appears_in_replicates)
