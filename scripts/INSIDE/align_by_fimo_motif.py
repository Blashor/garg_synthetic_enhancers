bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

tsv_lines = []
# INSIDE_256_g117.out
with open("INSIDE_jun13.txt") as file:
    # with open("INSIDE_klf4.txt") as file:
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
with open("fimo_40_160/fimo.tsv") as file:
    # with open("fimo_klf4/fimo.tsv") as file:
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
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = i_line + "_" + surv_rep
                # print(cell)
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                line_seq_holder[line_name][gen].append((motif, m_pos, cell[-1].strip()))
bag_of_motifs = ["SOX2", "NANOG", "ZIC3", "BACH2", "LEF1", "CTCF"]
for motif_of_choice in bag_of_motifs:
    os.system(f"mkdir -p weblogos/{motif_of_choice}")
    # os.system(f"mkdir -p weblogos/{motif_of_choice}_no_title")
    # os.system(f"mkdir -p weblogos/{motif_of_choice}_only_gens")
    dist_from_start = 32
    itr = 0
    evolved_motif_positions = {}
    for line in line_seq_holder:
        for gen in line_seq_holder[line]:
            # print(gen, line_seq_holder[line][gen][1:])

            if gen == "G150":
                score = line_seq_holder[line][gen][0][1]
                seq = line_seq_holder[line][gen][0][0]
                motifs = line_seq_holder[line][gen][1:]
                motifs.sort(key=lambda x: x[1])
                for motif1 in motifs:
                    if motif_of_choice == motif1[0]:
                        aligned_seq = seq[motif1[1] - dist_from_start : motif1[1] + dist_from_start]

                        if len(aligned_seq) == dist_from_start * 2:
                            itr += 1
                            # print(line, motif1[1], ,motif1[2], aligned_seq)
                            if line not in evolved_motif_positions:
                                evolved_motif_positions[line] = []
                            evolved_motif_positions[line].append(motif1[1])

    # This tracks proto positions and created motifs from the previous thing
    for gen in list(map(lambda x: f"G{x}", [1] + list(range(30, 310, 30)))):
        # for gen in list(map(lambda x: f"G{x}", range(5, 40, 5))):
        fasta_lines = []
        for line in evolved_motif_positions:
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            for itr, m_pos in enumerate(evolved_motif_positions[line]):
                aligned_seq = seq[m_pos - dist_from_start : m_pos + dist_from_start]
                fasta_lines.append(f">{line}_{itr}\n{aligned_seq}\n")
        with open("weblogo.fa", "w") as file:
            file.writelines(fasta_lines)
        gen_nice = gen.replace("G", "Generation:\\ ")
        os.system(
            f"weblogo --aspect-ratio 30 --color-scheme classic --stacks-per-line 100 --resolution 300 --format PDF --title {gen_nice} < weblogo.fa > weblogos/{motif_of_choice}/{motif_of_choice}_{gen}.pdf"
        )
        """
        os.system(
            f"weblogo --aspect-ratio 20 --color-scheme classic --stacks-per-line 100 --resolution 300 --format PNG < weblogo.fa > weblogos/{motif_of_choice}_no_title/{motif_of_choice}_{gen}.png"
        )
        os.system(
            f"weblogo --aspect-ratio 20 --color-scheme classic --stacks-per-line 100 --resolution 300 --format PNG --title {gen_nice} < weblogo.fa > weblogos/{motif_of_choice}_only_gens/{motif_of_choice}_{gen}.png"
        )
        """
