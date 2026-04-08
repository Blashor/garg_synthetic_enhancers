bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np


def smooth_array(data_array, window_size=4):
    smoothed = np.zeros_like(data_array, dtype=float)
    n = len(data_array)

    for i in range(n):
        # Calculate window bounds, making sure they stay within data_array
        start = max(0, i - window_size // 2)
        end = min(n, i + window_size // 2 + 1)
        smoothed[i] = np.mean(data_array[start:end])

    return smoothed


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
        "SP4",
    ]
)


allowed_motifs = set(
    [
        "PO5F1",
        "SOX2",
        "NANOG",
        "BACH2",
        "LEF1",
        "ZIC3",
        "CTCF",
    ]
)

palette_big = sns.color_palette("hls", 24)
palette_small = sns.color_palette("hls", 7)
motif_color_map = {
    "PO5F1": palette_big[0],
    "SOX2": palette_big[1],
    "NANOG": palette_big[2],
    "BACH2": palette_big[18],
    "LEF1": palette_big[16],
    "ZIC3": palette_small[2],
    "CTCF": palette_small[4],
}

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
            # print(gen)
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                # motif = f"{motif}_{cell[5]}"

                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                sign = cell[5]
                line_seq_holder[line_name][gen].append((motif, m_pos, sign))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import random

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import math

motif_colors = {}
palette = sns.color_palette("hls", 7)
color_idx = 0

fail_to_rise = set()
cached_data = {}

# --------- First Loop: Process and Cache Data ---------
for line in line_seq_holder:
    motif_set = set()
    score_arr = []
    m_gen = {}
    gen_list = sorted(line_seq_holder[line].keys(), key=lambda g: int(g.replace("G", "")))
    gen_to_idx = {}

    for i, gen in enumerate(gen_list):
        gen_int = int(gen.replace("G", ""))
        gen_to_idx[gen_int] = i

        score = float(line_seq_holder[line][gen][0][1])
        seq = line_seq_holder[line][gen][0][0]
        motifs = line_seq_holder[line][gen][1:]
        motifs.sort(key=lambda x: x[1])
        score = math.sqrt(score)
        score_arr.append(score)

        for m in motifs:
            m_id = f"{m[0]}{m[1]}{m[2]}"
            if m[0] in allowed_motifs:
                if m_id not in motif_set:
                    motif_set.add(m_id)
                    if gen_int not in m_gen:
                        m_gen[gen_int] = []
                    m_gen[gen_int].append(m[0])

                    if m[0] not in motif_colors:
                        motif_colors[m[0]] = palette[color_idx % len(palette)]
                        color_idx += 1

        if gen_int == 200:
            # > 147.5
            if score < 100 and len(motifs) > 0:
                fail_to_rise.add(line)

    # Cache only necessary data
    if line in fail_to_rise:
        cached_data[line] = {"score_arr": score_arr, "m_gen": m_gen, "gen_to_idx": gen_to_idx}

# --------- Second Loop: Plot with Subplots ---------
n = len(cached_data)
print(n)
fig, axs = plt.subplots(1, 1, figsize=(6, 4), sharex=True)


axs = [axs] * n  # Ensure iterable
import matplotlib.patches as patches

for ax, (line, data) in zip(axs, cached_data.items()):
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.tick_params(axis="both", width=2)
    score_arr = data["score_arr"]
    m_gen = data["m_gen"]
    gen_to_idx = data["gen_to_idx"]

    x = list(range(len(score_arr)))
    ax.plot(x, score_arr, label=f"{line}", color="grey", alpha=0.7, lw=2)

    # Arrows for motif appearance
    for gen_int in m_gen:
        idx = gen_to_idx.get(gen_int)
        if idx is None or idx >= len(score_arr):
            continue
        y = score_arr[idx]
        for motif in m_gen[gen_int]:
            color = motif_color_map.get(motif, "black")
            ax.annotate(
                "",
                xy=(idx, y + 0.5),
                xytext=(idx, y + 10),
                arrowprops=dict(
                    facecolor=color, shrink=0.00, width=0, headwidth=5, headlength=5, linewidth=0, alpha=0.9
                ),
            )
    # Legend
    legend_patches = [
        mpatches.FancyArrow(0, 0, 0.1, 0, color=motif_color_map[motif], label=motif.replace("O5", "OU5"))
        for motif in motif_colors
    ]
    ax.legend(handles=legend_patches, title="Motifs", bbox_to_anchor=(1.05, 1), loc="upper left")
    # ax.set_xlabel("Generation Index")
    # ax.set_ylabel("sqrt(Score)")
    # ax.set_title(f"Line: {line}")

plt.tight_layout()
plt.show()
