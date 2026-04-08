surv = 4
diff = False
line_seq_holder = {}
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


#
#
#

surv = 1


def generate_fimo():
    forty_fix = 0
    gen_add = 0
    # INSIDE_256_g117.out
    # with open("../data/INSIDE_jun13.txt") as file:
    with open(INSIDE_file) as file:
        current_stat = 0
        gen_tracker = {}
        replicates_mean = {}
        itr = 0
        first_lines = True
        for line in file:
            line = line.strip()
            if "NEW_FILE" in line:
                gen_add += 1

            if "" not in line:
                continue
            if "Gen:" not in line:
                continue
            cell = line.split("\t")
            cell[0] = "Gen: " + str(int(cell[0].split("Gen:")[-1]) + gen_add)
            current_stat += float(cell[3])
            if cell[1] not in line_seq_holder:
                line_seq_holder[cell[1]] = {}
            if cell[0] not in line_seq_holder[cell[1]]:
                line_seq_holder[cell[1]][cell[0]] = []
            if "', '" in cell[2]:
                cell[2] = cell[2].replace("(", "").replace(")", "").replace(" ", "").replace(",", "").replace("'", "")
            if "[" in cell[2]:
                print(cell[2])
            if itr % surv == 0:
                print(cell[1])
                line_seq_holder[cell[1]][cell[0]].append((cell[2], cell[3]))
            itr += 1

    fasta_seqs = []
    for line in line_seq_holder:
        for gen in line_seq_holder[line]:
            p_itr = 1
            for seq in line_seq_holder[line][gen]:
                # if gen == "Gen: 27":
                seq_id = line + "_" + gen.replace("Gen: ", "G") + "_" + str(p_itr)
                print(seq_id)
                p_itr += 1
                s_line = ">" + seq_id + "_" + seq[1] + "\n" + seq[0] + "\n"
                fasta_seqs.append(s_line)

    with open("INSIDE_meme.fasta", "w") as file:
        file.writelines(fasta_seqs)
    os.system(
        "export PATH=/opt/local/bin:/opt/local/libexec/meme-5.5.5:$PATH && fasta-get-markov -dna INSIDE_meme.fasta > markov_background.txt"
    )
    os.system(
        f"fimo --bfile markov_background.txt --no-pgc --oc {fimo_folder} ../data/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme INSIDE_meme.fasta"
    )
    """
    os.system(
        f"meme INSIDE_meme.fasta -dna -oc {fimo_folder} -nostatus -time 14400 -mod anr -nmotifs 6 -minw 5 -maxw 20 -objfun classic -revcomp -markov_order 0"
    )
    os.system(
        f"tomtom -verbosity 1 -oc {fimo_folder} -min-overlap 5 -dist ed -evalue -thresh 1 -no-ssc {fimo_folder}/meme.xml HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"
    )
    """


def motif_by_gen(gen_max=400):
    motif_gen_obj = {}
    motif_by_gen_lines = []
    with open(f"{fimo_folder}/fimo.tsv") as file:
        for line in file:
            cell = line.strip().split("\t")
            if len(cell) > 3:
                print(cell)
                if "line" in cell[2] or "pop" in cell[2] or "sequence" not in cell[2]:
                    # print(cell)
                    # line_rep, gen, p_num, coop_stat = cell[2].split("_")
                    # line, rep = line_rep.split("rep")
                    line, rep, gen, p_num, coop_stat = cell[2].split("_")

                    motif = cell[0].split("_")[0]
                    pvalue = float(cell[-3])
                    # print(line, gen[1:], motif)
                    if motif not in motif_gen_obj:
                        motif_gen_obj[motif] = [0] * gen_max
                    if pvalue < 10**-5:
                        motif_gen_obj[motif][int(gen[1:]) - 1] += 1
    for motif in motif_gen_obj:
        motif_by_gen_lines.append(motif + "\t" + "\t".join(map(str, motif_gen_obj[motif])) + "\n")
        # print(motif, "\t".join(map(str, motif_gen_obj[motif])), sep="\t")
    with open(f"{fimo_folder}/motif_by_gen.tsv", "w") as file:
        file.writelines(motif_by_gen_lines)


def meta_motif_by_gen():
    motif_gen_obj = {}
    os.system(f"mkdir -p {fimo_folder}/meta_motif")
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
                    if float(cell[7]) < pval and gen < 300:
                        if motif not in motif_gen_obj:
                            motif_gen_obj[motif] = {}
                        if gen not in motif_gen_obj[motif]:
                            motif_gen_obj[motif][gen] = {}
                        for coord in coords:
                            if coord not in motif_gen_obj[motif][gen]:
                                motif_gen_obj[motif][gen][coord] = [score]
                            else:
                                motif_gen_obj[motif][gen][coord].append(score)
    for motif in motif_gen_obj:
        for gen in motif_gen_obj[motif]:
            for coord in motif_gen_obj[motif][gen]:
                motif_gen_obj[motif][gen][coord] = len(motif_gen_obj[motif][gen][coord])
    for motif in motif_gen_obj:
        if motif in allowed_motifs:
            # print(motif)
            # print(motif_gen_obj[motif])
            row_indices = sorted(motif_gen_obj[motif].keys())
            col_indices = sorted(set(col for row in motif_gen_obj[motif].values() for col in row.keys()))
            # print(row_indices, col_indices)
            region_size = 512
            gen_num = max(row_indices)
            print(gen_num)
            matrix_data = np.zeros((gen_num, region_size))
            for gen in motif_gen_obj[motif]:
                for pos in motif_gen_obj[motif][gen]:
                    old_pos = pos
                    if pos > 511:
                        pos -= 512
                    matrix_data[gen - 1, pos - 1] += motif_gen_obj[motif][gen][old_pos]
            plt.figure()
            sns.heatmap(matrix_data, cmap="magma")
            plt.title(motif)
            plt.xlabel("Position (bp)")
            plt.ylabel("Generation")
            # print(matrix_data.shape)
            x_tick_positions = np.arange(128, 512, 128)  # Adjust interval as needed
            # x_tick_positions = np.array([50, 175, 300, 400, 525, 650])
            # x_tick_positions = np.array([0, 500, 1000])
            plt.xticks(ticks=x_tick_positions, labels=x_tick_positions, rotation=0, ha="center")
            y_tick_positions = [1] + list(np.arange(30, 300, 30))  # Adjust interval as needed
            plt.yticks(ticks=y_tick_positions, labels=y_tick_positions, rotation=0)
            # plt.axvline(x=matrix_data.shape[1] / 2, color="white", linewidth=2)
            plt.savefig(f"{fimo_folder}/meta_motif/" + motif + "_512.svg", dpi=300)
            # plt.show()


pval = 10**-6
allowed_motifs = set(["REST", "SOX2", "ZIC3", "NANOG", "OCT4", "CTCF"])
INSIDE_file = "INSIDE_1000.txt"
fimo_folder = "../data/INSIDE_40_definitive"
generate_fimo()
# motif_by_gen(36)
meta_motif_by_gen()
