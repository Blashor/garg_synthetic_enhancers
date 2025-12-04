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
            sns.lineplot(x=np.arange(len(plot_motif)), y=plot_motif, label=motif)
    plt.ylabel("Motif Occurences")
    plt.xlabel("Generation")
    plt.show()
    # print(motif_gen_obj[motif])


allowed_motifs = set(["BACH2", "SOX2", "ZIC3", "PO2F2", "PO5F1", "NANOG", "LEF1", "PO2F1", "CTCF", "REST"])
INSIDE_file = "INSIDE_jun13.txt"
fimo_folder = "/Users/blake/Documents/gargLab/INSIDE/INSIDE_40_definitive"
# generate_fimo()
# motif_by_gen(36)
meta_motif_by_gen()
