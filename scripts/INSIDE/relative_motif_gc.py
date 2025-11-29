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
with open("INSIDE_40_definitive/fimo.tsv") as file:
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
                line_seq_holder[line_name][gen].append([motif, m_pos, sign])
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

pseudo_motifs = {}
gen_skip = 1
for line in line_seq_holder:
    print(line)
    for gen in line_seq_holder[line]:
        gen_int = int(gen.replace("G", ""))
        gen_mod = gen_int % gen_skip
        # print(gen)
        if gen_int == 150:
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            motifs = line_seq_holder[line][gen][1:]
            motifs.sort(key=lambda x: x[1])
            # print(gen, motifs)
            pseudo_motifs[line] = motifs
            # for i1, motif1 in enumerate(motifs):
            #    print(line, motif1)


ecdf_holder = {}
gen_skip = 1
for line in line_seq_holder:
    for gen in line_seq_holder[line]:
        gen_int = int(gen.replace("G", ""))
        gen_mod = gen_int % gen_skip
        # print(gen)
        if gen_int == 1 or gen_int == 10 or gen_int == 50 or gen_int == 100:
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            motifs = line_seq_holder[line][gen][1:]
            motifs.sort(key=lambda x: x[1])
            # print(gen, motifs)
            for i1, motif1 in enumerate(pseudo_motifs[line]):
                # print(i1, motif1)
                if motif1[2] == "+":
                    if motif1[0] not in ecdf_holder:
                        ecdf_holder[motif1[0]] = {}
                    if gen_int not in ecdf_holder[motif1[0]]:
                        ecdf_holder[motif1[0]][gen_int] = {}
                        ecdf_holder[motif1[0]][gen_int] = {
                            "gc1": np.zeros([512]),
                            "total1": np.zeros([512]),
                            "gc2": np.zeros([512]),
                            "total2": np.zeros([512]),
                        }
                    if motif1[1] < 512:
                        seq = seq[:512]
                        first = seq[: motif1[1]]
                        second = seq[motif1[1] :]
                        # print(motif1[1], len(first), len(second))

                        # Add first (in reverse) to gc1 and total1
                        for i, letter in enumerate(reversed(first)):
                            if letter in "GC":
                                ecdf_holder[motif1[0]][gen_int]["gc1"][511 - i] += 1
                            ecdf_holder[motif1[0]][gen_int]["total1"][511 - i] += 1

                        # Add second (normal) to gc2 and total2
                        for i, letter in enumerate(second):
                            if letter in "GC":
                                ecdf_holder[motif1[0]][gen_int]["gc2"][i] += 1
                            ecdf_holder[motif1[0]][gen_int]["total2"][i] += 1
                    else:
                        seq = seq[512:]
                        mi = motif1[1] - 512
                        first = seq[:mi]
                        second = seq[mi:]
                        # print(motif1[1], len(first), len(second))

                        # Add first (in reverse) to gc1 and total1
                        for i, letter in enumerate(reversed(first)):
                            if letter in "GC":
                                ecdf_holder[motif1[0]][gen_int]["gc1"][511 - i] += 1
                            ecdf_holder[motif1[0]][gen_int]["total1"][511 - i] += 1

                        # Add second (normal) to gc2 and total2
                        for i, letter in enumerate(second):
                            if letter in "GC":
                                ecdf_holder[motif1[0]][gen_int]["gc2"][i] += 1
                            ecdf_holder[motif1[0]][gen_int]["total2"][i] += 1
os.system("mkdir -p relative_motif_gc_within_across_INSIDE_def")
# print(ecdf_holder)
# Gen 1-30, 31-60, 61-90, 91-120,121-150
for motif_id in ecdf_holder:
    print(motif_id[:-2])
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.tick_params(axis="both", width=2)

    for gen_int in ecdf_holder[motif_id]:
        gc1 = ecdf_holder[motif_id][gen_int]["gc1"]
        gc2 = ecdf_holder[motif_id][gen_int]["gc2"]
        total1 = ecdf_holder[motif_id][gen_int]["total1"]
        total2 = ecdf_holder[motif_id][gen_int]["total2"]
        print(total1)
        # Concatenate the arrays
        gc = np.concatenate([gc1, gc2])  # shape (1024,)
        total = np.concatenate([total1, total2])

        gc10 = ecdf_holder[motif_id][1]["gc1"]
        gc20 = ecdf_holder[motif_id][1]["gc2"]
        total10 = ecdf_holder[motif_id][1]["total1"]
        total20 = ecdf_holder[motif_id][1]["total2"]
        gc0 = np.concatenate([gc10, gc20])  # shape (1024,)
        total0 = np.concatenate([total10, total20])
        print(np.sum(gc))
        if np.sum(gc) > 0:
            # Avoid division by zero
            total0[total0 == 0] = 1
            total[total == 0] = 1
            # Compute frequency (or normalized values)
            data_array = smooth_array(gc / total, 10) - smooth_array(gc0 / total0, 10)  # shape (1024,)
            print(data_array)
            data_array = data_array[256 : 512 + 256]
            # Now plot the histogram

            # Plot histogram
            x_values = list(range(-256, -256 + len(data_array)))
            plt.plot(
                x_values, data_array, color="black", alpha=(gen_int / 167) + 0.1, linewidth=1, label=f"Gen: {gen_int}"
            )
            plt.ylim([-0.1, 0.2])
            plt.xlabel("Relative Distance (bp)")
            plt.ylabel("GC Change")
            plt.tight_layout()
            plt.title(motif_id)
            plt.xticks(range(-256, 257, 128))
    plt.legend()
    plt.show()
    plt.savefig(f"relative_motif_gc_within_across_INSIDE_def/{motif_id}.svg", dpi=300)
    plt.close()
