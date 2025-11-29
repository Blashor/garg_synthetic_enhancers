bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os


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
            print(gen, motifs)
            for i1, motif1 in enumerate(motifs):
                for i2, motif2 in enumerate(motifs):
                    if (
                        motif1[0].split("_")[0] in allowed_motifs
                        and motif2[0].split("_")[0] in allowed_motifs
                        and i1 != i2
                    ):
                        motif_id = f"{motif1[0]}-{motif2[0]}"
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
