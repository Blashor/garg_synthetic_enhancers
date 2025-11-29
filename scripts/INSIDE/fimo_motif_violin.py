bed_lines = []
gen_add = 0
survivors = 1
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def set_style(figsize=(10, 6)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


set_style()
tsv_lines = []
# INSIDE_256_g117.out
INSIDE_file = "INSIDE_1000.txt"  # "INSIDE_jun13.txt"
fimo_folder = "/Users/blake/Documents/gargLab/INSIDE/INSIDE_1000"  # fimo_40_160
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
            # print(cell[0])
            line_seq_holder[cell[1]][cell[0]] = []
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
# print(line_seq_holder)
with open(f"{fimo_folder}/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            # print(cell[2])
            g_line, i_line, oth_rep, gen, surv_rep, score = cell[2].split("_")
            # i_line, gen, surv_rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            gen_n = int(gen.split("G")[-1])
            # print(gen)
            # print(cell)
            if float(cell[-3]) < 10**-6:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = g_line + "_" + i_line + "_" + oth_rep + f"_{surv_rep}"
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                # print(line_seq_holder[line_name])
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

violin_dict = {}
color_motif = {}
clrs = sns.color_palette("deep", 10) + sns.color_palette("husl", 191)
clr_i = 0
for line in line_seq_holder:
    for gen in line_seq_holder[line]:
        gen_num = int(gen.replace("G", ""))
        last_gen = f"G{gen_num-1}"
        if gen != "G1" and gen_num < 300:
            try:
                score = math.sqrt(float(line_seq_holder[line][gen][0][1]))
                last_score = math.sqrt(float(line_seq_holder[line][last_gen][0][1]))
                seq = line_seq_holder[line][gen][0][0]
                motifs = line_seq_holder[line][gen][1:]
                print(score, motifs)
                motif_holder = {}
                for motif, m_pos in motifs:
                    if motif not in motif_holder:
                        motif_holder[motif] = 0
                    motif_holder[motif] += 1
                for motif in motif_holder:
                    motif_id = f"{motif}_{motif_holder[motif]}"
                    if motif not in color_motif:
                        color_motif[motif] = clrs[clr_i]
                        clr_i += 1
                    if motif_id not in violin_dict:
                        violin_dict[motif_id] = []
                    violin_dict[motif_id].append(score)
            except Exception:
                pass
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


for violin in violin_dict:
    if np.median(violin_dict[violin]) > 0 and len(violin_dict[violin]) > 600 and violin.split("_")[0] in allowed_motifs:
        # print(violin, len(violin_dict[violin]))
        motif, m_num = violin.split("_")
        clr = color_motif[motif]
        box_clusters.append((f"{motif} ({m_num})", violin_dict[violin], np.median(violin_dict[violin]), clr))

box_clusters.sort(key=lambda x: x[2], reverse=True)
labels, violins, med, clr = zip(*box_clusters)

ax = sns.boxplot(data=violins, width=0.5, linewidth=0.5, fliersize=0)
ax.set_xticklabels(labels, rotation=45, ha="right")
for patch, color in zip(ax.patches, clr):
    patch.set_facecolor(color)
plt.show()
# with open("fimo_40/INSIDE.tsv", "w") as file:
#    file.writelines(tsv_lines)
