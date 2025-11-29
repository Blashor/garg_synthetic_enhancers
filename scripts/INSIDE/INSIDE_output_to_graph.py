import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy
import numpy as np

labels = {
    "line1": "SNP 0.4%",
    "line2": "SNP 1%",
    "line3": "1-3 Duplications",
    "line4": "1-3 Duplications + SNP 1%",
    "line5": "SNP 5%",
    "line6": "SNP 20%",
}
# sns.set_style("whitegrid")


def main():
    # "INSIDE_16693679.out"
    "INSIDE_16879235.out"
    "INSIDE_2reps.out"
    "gen50_1.out"
    # INSIDE_10len.out
    # "INSIDE_mar.out"
    # "INSIDE_22b.out"
    #
    gen_add = 0
    surv = 4
    diff = False
    file_n = "/Users/blake/Downloads/INSIDE_genomeC.out"
    # file_n = "INSIDE_graded_on4_3model.txt"
    # file_n = "INSIDE_jun17.txt"
    # file_n = "INSIDE_klf4.txt"
    # file_n = "INSIDE_33786447.out"
    # file_n = "INSIDE_genome_sep12.txt"
    # file_n = "INSIDE_shrinking_window.txt"
    with open(file_n) as file:
        current_stat = 0
        gen_tracker = {}
        replicates_mean = {}
        itr = 0
        first_lines = True
        for line in file:
            # print(line)
            # mut
            # center
            # _
            if "NEW_FILE" in line:
                gen_add += 1
            if "Gen:" not in line:
                continue
            # if "populationATAC" not in line:
            #    continue
            if "line256rep6" in line:
                print(line)
            cell = line.split("\t")
            current_stat = float(cell[3])
            # 4, 3
            cell[1] = cell[1].replace("population", "")
            if itr % surv == surv - 1:
                if cell[1] not in gen_tracker:
                    gen_tracker[cell[1]] = []
                gen_tracker[cell[1]].append(current_stat)
                current_stat = 0

            itr += 1

        for line in gen_tracker:
            print(np.sqrt(gen_tracker[line]))
            if line[:-1] not in replicates_mean:
                replicates_mean[line[:-1]] = gen_tracker[line]
            else:
                for itr2 in range(len(replicates_mean[line[:-1]])):
                    replicates_mean[line[:-1]][itr2] += gen_tracker[line][itr2]
            gen_tracker[line] = gen_tracker[line][:10000]

            sns.lineplot(
                x=range(len(gen_tracker[line])),
                y=np.sqrt(gen_tracker[line]),
                label=f"{line}_{round(np.sqrt(gen_tracker[line][0]),1)}->{round(np.sqrt(gen_tracker[line][-1]),1)}",
            )
            print(line, gen_tracker[line][-1])
        print(replicates_mean)
        for line in replicates_mean:
            for itr2 in range(len(replicates_mean[line])):
                replicates_mean[line][itr2] = replicates_mean[line][itr2]
            if diff == True:
                # replicates_mean[line] = replicates_mean[line][:300]
                ax = sns.lineplot(
                    x=range(len(numpy.diff(replicates_mean[line]))),
                    y=numpy.diff(replicates_mean[line]),
                    label=line,
                )
            else:
                print(line)

                # line_l = line.replace("line_margin_", "Margin Length: ")
                """
                sns.lineplot(
                    x=range(len(replicates_mean[line])),
                    y=replicates_mean[line],
                    label=line,
                )
                """

        # plt.ylabel("Odds Ratio")
        # plt.xlabel("Generation #")

        # plt.legend(loc="upper left")
        # plt.title(file_n)
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(f"{file_n}.png", dpi=600)
        plt.show()


main()
"""
import json

line_obj = {}
line_arr = []
with open("INSIDE_16693679.out") as file:
    for line in file:
        # print(line)
        cell = line.split("\t")
        if cell[0] == "Gen: 43":
            if cell[1] not in line_obj:
                line_obj[cell[1]] = []
            line_obj[cell[1]].append((cell[2], float(cell[3])))

for lines in line_obj:
    line_r = line_obj[lines]
    line_arr.append(line_r)

with open("latest.checkpoint", "w") as file:
    file.write(json.dumps({"scored_parents": line_arr, "current_gen": 43}))

with open("latest.checkpoint") as file:
    obj = json.load(file)
    print(obj)
"""
