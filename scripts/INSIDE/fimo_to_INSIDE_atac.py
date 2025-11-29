bed_lines = {}
offset = 2448  # enhancer insertion site
with open("fimo_40/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, rep, score = cell[2].split("_")
            if i_line == "line256rep6":
                # chr_n = ord(i_line[-1]) - ord("a") + 1
                chr_n = 1
                chrom = f"chr{chr_n}"
                # print(gen)
                if gen not in bed_lines:
                    bed_lines[gen] = []
                    # print(cell)
                if rep == "1":
                    if float(cell[-3]) < 10**-5:
                        print(cell[0].split("_")[0])
                        cell[0] = cell[0].split("_")[0]  # + "_" + score
                        b_line = (
                            "\t".join(
                                [chrom, str(int(cell[3]) + offset), str(int(cell[4]) + offset), cell[0], "1", cell[5]]
                            )
                            + "\n"
                        )
                        print(b_line)
                        bed_lines[gen].append(b_line)
# print(bed_lines)
for gen in ["G1", "G30", "G60", "G90", "G120", "G150"]:
    with open(f"fimo_40/INSIDE_atac_fimo_{gen}.bed", "w") as file:
        file.writelines(bed_lines[gen])
