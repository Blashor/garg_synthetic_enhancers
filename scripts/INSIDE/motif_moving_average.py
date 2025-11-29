bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

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
with open("INSIDE_jun17.txt") as file:
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
                line_seq_holder[line_name][gen].append((motif, m_pos, cell[5]))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

motif_ma = {}
gens_to_meta = ["G200"]
ecdf_holder = {}
gen_skip = 1
for gen in gens_to_meta:
    motif_by_pos = np.zeros(512)
    motif_by_pos_zic = np.zeros(512)
    for s_i, sample in enumerate(line_seq_holder):
        score = line_seq_holder[sample][gen][0][1]
        seq = line_seq_holder[sample][gen][0][0]
        motifs = line_seq_holder[sample][gen][1:]
        motifs.sort(key=lambda x: x[1])
        if s_i < 10000:
            for m in motifs:
                name, pos, sense = m
                if "NANOG" in name or "PO5F1" in name or "SOX2" in name:
                    if pos >= 512:
                        motif_by_pos[pos - 512] += pos
                    else:
                        motif_by_pos[pos] += pos
                if "ZIC3" in name:
                    if pos >= 512:
                        motif_by_pos_zic[pos - 512] += pos
                    else:
                        motif_by_pos_zic[pos] += pos
    window = 8
    motif_ma = np.convolve(motif_by_pos, np.ones(window) / window, mode="valid")

    motif_ma2 = np.convolve(motif_by_pos_zic, np.ones(window) / window, mode="valid")
    motif_ma = motif_ma / max(motif_ma)
    motif_ma2 = motif_ma2 / max(motif_ma2)
    sns.lineplot(motif_ma, label=f"{gen}_OSN", alpha=0.8, linewidth=1.25)
    sns.lineplot(motif_ma2, label=f"{gen}_ZIC3", alpha=0.8, linewidth=1.25, linestyle="--")
plt.axvline(x=128, color="b", linestyle="--")
plt.axvline(x=512 - 128, color="b", linestyle="--")
plt.axhline(y=0.5, color="black", linestyle="--")
plt.title("ZIC")
plt.show()


# print(line_seq_holder)
