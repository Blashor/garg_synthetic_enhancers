bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# sns.set_style("whitegrid")
tsv_lines = []
# INSIDE_256_g117.out

INSIDE_file = "../data/INSIDE_jun13.txt"
INSIDE_file = "../data/INSIDE_graded_on4_3model.txt"
fimo_folder = "../data/INSIDE_40_definitive"  # fimo_40_160 "INSIDE_40_definitive/fimo.tsv"
# INSIDE_file = "INSIDE_1000.txt"
# INSIDE_file = "../data/INSIDE_graded_on4_3model.txt"
# fimo_folder = "../data/INSIDE_1000"  # fimo_40_160 "INSIDE_40_definitive/fimo.tsv"
with open(INSIDE_file) as file:
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
            # print(cell[1], cell[0])
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
# print(line_seq_holder)
with open(f"{fimo_folder}/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            # g_line, i_line, oth_rep, gen, surv_rep, score = cell[2].split("_")
            i_line, gen, surv_rep, score = cell[2].split("_")
            chr_n = ord(i_line[-1]) - ord("a") + 1
            # chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            gen_n = int(gen.split("G")[-1])
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                # line_name = g_line + "_" + i_line + "_" + oth_rep + f"_{surv_rep}"
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                # print(line_name, gen, line_seq_holder[line_name])
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

violin_dict = {}
motif_score_dict = {}
color_motif = {}
clrs = sns.color_palette("deep", 10) + sns.color_palette("husl", 300)
clr_i = 0
replicate_num_motif_occurs = {}
score_dict = {}
for line in line_seq_holder:
    sample, surv_num = line.split("_")
    for gen in line_seq_holder[line]:
        gen_num = int(gen.replace("G", ""))
        gen_skip = 1
        last_gen = f"G{gen_num-gen_skip}"
        if sample not in score_dict:
            score_dict[sample] = {}
        if gen not in score_dict[sample]:
            score_dict[sample][gen] = []

        if gen_num > gen_skip and gen_num % gen_skip == 0:
            # print(gen, last_gen)
            # try:
            score = math.sqrt(float(line_seq_holder[line][gen][0][1]))
            score_dict[sample][gen].append(score)
            last_score = math.sqrt(float(line_seq_holder[line][last_gen][0][1]))
            score_diff = gen_num
            seq = line_seq_holder[line][gen][0][0]
            motifs = set(line_seq_holder[line][gen][1:])
            # print(motifs)
            if len(motifs) == 0:
                pass
                # motifs.add(("None", 0))
            last_motifs = set(line_seq_holder[line][last_gen][1:])
            new_motifs = motifs - last_motifs
            good_new_motifs = set()
            for motif in new_motifs:
                motif_name = motif[0]
                good_new_motifs.add(motif_name)
            motif_holder = {}

            for motif, m_pos in motifs:
                if motif in good_new_motifs:
                    if motif not in motif_holder:
                        motif_holder[motif] = 0
                    motif_holder[motif] += 1
            for motif in motif_holder:
                motif_id = f"{motif}_{motif_holder[motif]}"
                # print(motif_id)
                if motif_id not in replicate_num_motif_occurs:
                    replicate_num_motif_occurs[motif_id] = set()
                if line not in replicate_num_motif_occurs[motif_id]:
                    replicate_num_motif_occurs[motif_id].add(line)
                if motif not in color_motif:
                    color_motif[motif] = clrs[clr_i]
                    clr_i += 1
                if motif_id not in violin_dict:
                    violin_dict[motif_id] = []
                violin_dict[motif_id].append(score_diff)
                if motif_id not in motif_score_dict:
                    motif_score_dict[motif_id] = []
                motif_score_dict[motif_id].append((score))
        # except Exception:
        #    pass
        # holder = [line, gen, score, seq, "\t".join(map(lambda x: str(x[0]) + "\t" + str(x[1]), motifs))]
        # tsv_lines.append("\t".join(holder) + "\n")
    # print(violin_dict)
box_clusters = []
# print(clr_i)
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

allowed_motifs = set(["BACH2", "SOX2", "ZIC3", "PO2F2", "PO5F1", "NANOG", "LEF1", "PO3F1", "PO2F1", "CTCF", "None"])

display_score = {}
for sample in score_dict:
    for i, gen in enumerate(score_dict[sample]):
        if sample not in display_score:
            display_score[sample] = []
            #

            # display_score[sample].append(np.mean(score_dict[sample][gen]))
        if i > 0:
            """
            score_met = (np.mean(score_dict[sample][gen]) - np.mean(score_dict[sample][f"G{i}"])) - (
                np.mean(score_dict[sample][f"G{i}"]) - np.mean(score_dict[sample][f"G{i-1}"])
            )
            """
            score_met = np.mean(score_dict[sample][gen])
            display_score[sample].append(score_met)


import matplotlib.colors as mcolors

cmap = sns.color_palette("Spectral_r", as_cmap=True)
norm = mcolors.Normalize(vmin=0, vmax=40)
norm_score = mcolors.Normalize(vmin=0, vmax=180)
for violin in violin_dict:
    # print(violin)
    if (np.median(violin_dict[violin]) > 0 and (len(replicate_num_motif_occurs[violin]) / 4) > 5) and violin.split("_")[
        0
    ] in allowed_motifs:
        if int(violin.split("_")[-1]) < 6:
            # print(violin, len(violin_dict[violin]))
            print(violin, len(replicate_num_motif_occurs[violin]) / 4)
            motif = violin.split("_")[0]
            motif_name = violin.split("_")
            motif_name[0] = motif_name[0].replace("PO", "POU")
            if motif_name[0] == "PO5F1":
                motif_name[0] = "POU5F1"
            if motif_name[0] == "PO2F1":
                motif_name[0] = "POU2F1"
            if motif_name[0] == "PO2F2":
                motif_name[0] = "POU2F2"
            motif_name = f"{motif_name[0]} ({motif_name[1]})"
            # clr = color_motif[motif]
            clr = cmap(norm(len(replicate_num_motif_occurs[violin]) / 4))
            clr_med = cmap(norm_score(np.median(motif_score_dict[violin])))
            box_clusters.append((motif_name, violin_dict[violin], np.median(violin_dict[violin]), clr, clr_med))

box_clusters.sort(key=lambda x: x[2], reverse=True)
labels, violins, med, clr, clr_med = zip(*box_clusters)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7), sharex=True, gridspec_kw={"height_ratios": [4, 7]})
# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax2 = ax1.twinx()
# Plot the line plot
for sample in display_score:
    print(sample)
    # window = 2
    # display_score[sample] = np.convolve(display_score[sample], np.ones(window) / window, mode="valid")
    color = "black"
    alpha = 0.5
    linewidth = 1
    if sample == "line256rep6":
        continue
    elif sample == "line256rep28":
        color = "#9cd7a4"
        linewidth = 2
        alpha = 1
    elif sample == "line256rep3":
        continue

    ax1.plot(
        range(len(display_score[sample])),
        display_score[sample],
        color=color,
        label=sample,
        alpha=alpha,
        linewidth=linewidth,
    )
for sample in display_score:
    print(sample)
    # window = 2
    # display_score[sample] = np.convolve(display_score[sample], np.ones(window) / window, mode="valid")
    color = "black"
    alpha = 0
    linewidth = 2
    if sample == "line256rep6":
        color = "#c03830"
        alpha = 1

    elif sample == "line256rep28":
        continue
    elif sample == "line256rep3":
        color = "#7192BE"
        alpha = 1
    else:
        continue

    ax1.plot(
        range(len(display_score[sample])),
        display_score[sample],
        color=color,
        label=sample,
        alpha=alpha,
        linewidth=linewidth,
    )

# Plot the boxplot horizontally
ax_o = sns.boxplot(data=violins, orient="h", width=0.5, linewidth=0.5, ax=ax2, fliersize=0)
# ax_o = sns.stripplot(data=violins, orient="h", ax=ax2, palette=clr, size=2)
"""
ax_o = sns.violinplot(
    data=violins,
    orient="h",
    bw_method="silverman",
    width=1,
    linewidth=0.5,
    ax=ax2,
    inner="point",
    density_norm="count",
    palette=clr,
)
"""
# Set colors for each box
for patch, color in zip(ax2.patches, clr_med):
    patch.set_facecolor(color)

# Set y-tick labels
ax2.set_yticklabels(labels, fontsize=7)
plt.xlabel("Generation")
plt.xlim(-1, len(display_score[sample]) + 1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_score)
sm.set_array([])
cbar = fig.colorbar(sm, orientation="horizontal", pad=0.12, ax=ax2, shrink=0.5)


def set_style(ax):
    # fig = plt.figure(figsize=figsize)
    # ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


set_style(ax2)
set_style(ax1)
fig.subplots_adjust(bottom=0.25)
# Adjust layout to prevent overlap
fig.tight_layout()
plt.savefig("Violin_lineplot.svg", dpi=200)
plt.show()
