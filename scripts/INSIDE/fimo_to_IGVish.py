bed_lines = {}
with open("fimo_40/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            i_line, gen, rep, score = cell[2].split("_")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            # print(gen)
            if gen not in bed_lines:
                bed_lines[gen] = []
                # print(cell)
            if rep == "1":
                if float(cell[-3]) < 10**-6:
                    print(cell[0].split("_")[0])
                    cell[0] = cell[0].split("_")[0]  # + "_" + score
                    b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                    print(b_line)
                    bed_lines[gen].append(b_line)
# print(bed_lines)
for gen in ["G1", "G30", "G60", "G90", "G120", "G150"]:
    with open(f"fimo_40/IGVish_{gen}.bed", "w") as file:
        file.writelines(bed_lines[gen])
