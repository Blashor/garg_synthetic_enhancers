import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# fimo_motif_line
def meta_motif_by_gen():
    motif_gen_obj = {}
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)
    with open(f"{fimo_folder}/fimo.tsv") as file:
        for line in file:
            cell = line.strip().split("\t")
            if len(cell) > 3:
                if "line" in cell[2] or "population" in cell[2]:
                    i_line, gen, p_num, coop_stat = cell[2].split("_")
                    gen = int(gen[1:])
                    motif = cell[0].split("_")[0]
                    score = float(cell[6])
                    p_val = float(cell[7])
                    coords = range(int(cell[3]), int(cell[4]))
                    if float(cell[7]) < 10**-6 and gen < 300:
                        if motif not in motif_gen_obj:
                            motif_gen_obj[motif] = {}
                        if gen not in motif_gen_obj[motif]:
                            motif_gen_obj[motif][gen] = []
                        motif_gen_obj[motif][gen].append(score)

    for motif in motif_gen_obj:
        plot_motif = np.zeros([300])
        if motif in allowed_motifs:
            for gen in motif_gen_obj[motif]:
                if gen < len(plot_motif):
                    plot_motif[gen] = len(motif_gen_obj[motif][gen]) / 40
            motif = motif.replace("PO", "POU")
            if motif in motif_color_map:
                sns.lineplot(x=np.arange(len(plot_motif)), y=plot_motif, label=motif, color=motif_color_map[motif])
            else:
                sns.lineplot(x=np.arange(len(plot_motif)), y=plot_motif, label=motif)
    plt.ylabel("Motif Occurences")
    plt.xlabel("Generation")
    plt.show()
    # print(motif_gen_obj[motif])


palette_big = sns.color_palette("hls", 24)
palette_small = sns.color_palette("hls", 7)
palette_massive = sns.color_palette("hls", 48)
motif_color_map = {
    "POU5F1": palette_big[0],
    "SOX2": palette_big[1],
    "NANOG": palette_big[2],
    "BACH2": palette_big[18],
    "LEF1": palette_big[16],
    "ZIC3": palette_small[2],
    "CTCF": palette_small[4],
    "POU2F2": palette_big[11],
    "POU3F1": palette_big[12],
    "KLF4": palette_massive[42],
    "SP2": palette_massive[40],
    "SP3": palette_massive[41],
    "None": (0, 0, 0, 0.9),
}


allowed_motifs = set(
    ["BACH2", "SOX2", "ZIC3", "PO2F2", "PO5F1", "NANOG", "LEF1", "PO2F1", "CTCF", "REST", "KLF4", "SP2", "SP3"]
)

INSIDE_file = "../data/data/INSIDE_jun13.txt"
fimo_folder = "../data/data/INSIDE_40_definitive"
# generate_fimo()
# motif_by_gen(36)
meta_motif_by_gen()
