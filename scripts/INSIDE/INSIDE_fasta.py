""

import os
import numpy as np

dir20 = "/Users/blake/Downloads/INSIDE20b"

gen_holder = {}

with open(f"INSIDE_jun13.txt") as file:
    for line in file:
        if "Gen:" in line:
            cell = line.split("\t")
            gen = int(cell[0].split("Gen: ")[-1]) - 1
            gene_line = cell[1]
            if gene_line not in gen_holder:
                gen_holder[gene_line] = {}
            if gen not in gen_holder[gene_line]:
                gen_holder[gene_line][gen] = cell
i = 0
lines = []
f_lines = []
for gene_line in gen_holder:
    # print(gene_line)
    for num in range(150):
        cell = gen_holder[gene_line][num]
        # new_lines = gen_holder[gene_line][num].replace("('", "").replace("')", "").split("', '")
        # gen, name_id, seq1 = new_lines[0].split("\t")
        # seq2, eor = new_lines[1].split("\t")
        f_lines.append(f">{cell[1]}_G{num+1}_S1\n{cell[2]}\n")
        # f_lines.append(f">{name_id}_G{num+1}_S2\n{seq2}\n")
        # new_line = "".join(new_lines)
        # lines.append(new_line + "\n")
    if max(gen_holder[gene_line].keys()) < 50:
        l20 = int(int(gene_line.split("_")[-1]) / 50)
        # print(i, l20, gene_line, max(gen_holder[gene_line].keys()))

        i += 1
# with open("INSIDE_1000.txt", "w") as file:
#    file.writelines(lines)

with open("INSIDE_40.fa", "w") as file:
    file.writelines(f_lines)
# 3, 4, 7,9, 14
