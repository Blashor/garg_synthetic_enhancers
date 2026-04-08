bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy.optimize import curve_fit

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
                motif = f"{motif}_{cell[5]}"

                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                line_seq_holder[line_name][gen].append((motif, m_pos))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

ecdf_holder = {}
gen_skip = 1
for line in line_seq_holder:
    for gen in line_seq_holder[line]:
        gen_mod = int(gen.replace("G", "")) % gen_skip
        # print(gen)
        if gen != "G30dhjfhhbf0":
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            motifs = line_seq_holder[line][gen][1:]
            motifs.sort(key=lambda x: x[1])
            for i1, motif1 in enumerate(motifs):
                # print(motif1)
                for i2, motif2 in enumerate(motifs):
                    if (
                        motif1[0].split("_")[0] in allowed_motifs
                        and motif2[0].split("_")[0] in allowed_motifs
                        and i1 != i2
                    ):
                        motif_id = f"{motif1[0].split("_")[0]}-{motif2[0].split("_")[0]}"
                        if motif_id not in ecdf_holder:
                            ecdf_holder[motif_id] = []
                        if motif1[1] > 512 and motif2[1] > 512:
                            pass
                        elif motif1[1] < 512 and motif2[1] < 512:
                            pass
                        else:
                            continue
                        rel_pos = motif1[1] - motif2[1]
                        ecdf_holder[motif_id].append(rel_pos)
#os.system("mkdir -p relative_sign_motif_within_across_INSIDE_def")
# print(ecdf_holder)

period_output = []
for motif_id in ecdf_holder:
    print(ecdf_holder[motif_id])
    if len(ecdf_holder[motif_id]) > 20 * 400:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
        ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
        ax.tick_params(axis="both", width=2)
        # Plot the histogram in the first subplot
        sns.histplot(data=ecdf_holder[motif_id], bins=int(512 / 8), alpha=1, color="black")

        # --- compute histogram for fitting ---
        counts, bins = np.histogram(ecdf_holder[motif_id], bins=int(512 / 8))
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # sine model
        amp = np.mean(counts)

        def sine_func(x, A, B, D):
            return -amp * np.cos(B * x) + D

        # initial guess (helps convergence)
        p0 = [np.max(counts), 2 * np.pi / 160, np.mean(counts)]

        # fit
        popt, _ = curve_fit(sine_func, bin_centers, counts, p0=p0)

        # smooth curve for plotting
        x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 900)
        y_fit = sine_func(x_fit, *popt)

        A, B, D = popt

        # Compute periodicity
        period = 2 * np.pi / B
        print(f"{motif_id} periodicity: {period:.2f} bp")
        period_output.append(f"{motif_id}\t{period:.2f}\n")
        # overlay fit
        #plt.plot(x_fit, y_fit, color="grey", linewidth=2)

        plt.xlabel("Relative Distance (bp)")
        plt.ylabel("Frequency")
        plt.title(motif_id + "_" + str(int(period)))
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"../figures/relative_motif_within_across_INSIDE_def/{motif_id}noplot.svg", dpi=300)
        #plt.close()
with open("../data/5c_relative_motif_sinefit.tsv","w") as file:
    file.writelines(period_output)
