bed_lines = []
gen_add = 0
survivors = 4
line_seq_holder = {}
import math

tsv_lines = []
# INSIDE_256_g117.out

# INSIDE_file = "INSIDE_jun13_160gens.txt"
# fimo_folder = "fimo_40_160"
fimo_file = "adipose_smooth_plus"
with open("/Users/blake/Downloads/output_plus.txt") as file:
    # with open("INSIDE_genome_sep12_75gens.txt") as file:
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
        # print(cell[1], cell[0])
        print(cell[0], cell[1])
        line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
# print(line_seq_holder)
with open(f"{fimo_file}/fimo.tsv") as file:
    for line in file:
        cell = line.split("\t")
        if cell[0] != "motif_id" and len(cell) > 5:
            if float(cell[-3]) < 10**-5:
                print(cell[0].split("_")[0], cell)
            # RED cell[2] = cell[2].replace("_rep", "rep")
            # print(cell[2])
            # print(cell[2])
            p, pop, i_line, gen, surv_rep, score = cell[2].split("_")
            # i_line, gen, surv_rep, score = cell[2].split("_")
            # RED i_line = i_line.replace("rep", "_rep")
            # chr_n = ord(i_line[-1]) - ord("a") + 1
            chr_n = i_line.split("rep")[-1]
            chrom = f"chr{chr_n}"
            # print(gen)
            # print(cell)
            if float(cell[-3]) < 10**-5:
                # print(cell[0].split("_")[0])
                cell[0] = cell[0].split("_")[0] + "_" + score
                motif, score = cell[0].split("_")
                b_line = "\t".join([chrom, cell[3], cell[4], cell[0], "1", cell[5]]) + "\n"
                line_name = f"{p}_{pop}_{i_line}_{surv_rep}"
                # print(line_name)
                m_pos = int((int(cell[3]) + int(cell[4])) / 2)
                # -math.log10(float(cell[-3]))
                line_seq_holder[line_name][gen].append((motif, int(cell[3]), int(cell[4])))
                # print(line_seq_holder[line_name][gen])
                # bed_lines.append(b_line)

gen_skip = 10
for line in line_seq_holder:
    for gen in line_seq_holder[line]:
        gen_mod = int(gen.replace("G", "")) % gen_skip
        if gen == "G1" or gen_mod == 0 or gen == "G300":
            score = line_seq_holder[line][gen][0][1]
            seq = line_seq_holder[line][gen][0][0]
            motifs = line_seq_holder[line][gen][1:]
            motifs.sort(key=lambda x: x[1])
            holder = [
                line,
                gen,
                score,
                seq,
                "\t".join(map(lambda x: str(x[0]) + "\t" + str(x[1]) + "\t" + str(round(x[2])), motifs)),
            ]
            # print(holder)
            tsv_lines.append("\t".join(holder) + "\n")

with open(f"{fimo_file}/INSIDE_motifs.tsv", "w") as file:
    file.writelines(tsv_lines)
