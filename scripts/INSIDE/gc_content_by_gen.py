import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy
import Bio.SeqUtils


labels = {
    "line1": "SNP 0.4%",
    "line2": "SNP 1%",
    "line3": "1-3 Duplications",
    "line4": "1-3 Duplications + SNP 1%",
    "line5": "SNP 5%",
    "line6": "SNP 20%",
}
sns.set_style("whitegrid")


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
    file_n = "INSIDE_jun13.txt"
    file_n = "INSIDE_graded_on4_3model.txt"
    file_n = "INSIDE_jun17.txt"
    # file_n = "INSIDE_33786447.out"
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
            if "" not in line:
                continue
            if "NEW_FILE" in line:
                gen_add += 1
            if "Gen:" not in line:
                continue
            if "line256rep29" in line:
                print(line)
            cell = line.split("\t")
            # print(cell[2])
            current_stat += float(cell[3])
            # 4, 3
            if itr % surv == surv - 1:
                if cell[1] not in gen_tracker:
                    gen_tracker[cell[1]] = [0]
                gen_tracker[cell[1]].append(Bio.SeqUtils.gc_fraction(cell[2]))
                current_stat = 0

            itr += 1

        for line in gen_tracker:
            print(gen_tracker[line])
            if line[:-1] not in replicates_mean:
                replicates_mean[line[:-1]] = gen_tracker[line]
            else:
                for itr2 in range(len(replicates_mean[line[:-1]])):
                    replicates_mean[line[:-1]][itr2] += gen_tracker[line][itr2]

            sns.lineplot(
                x=range(len(gen_tracker[line])),
                y=gen_tracker[line],
                # label=line,
            )

        print(replicates_mean)
        for line in replicates_mean:
            for itr2 in range(len(replicates_mean[line])):
                replicates_mean[line][itr2] = replicates_mean[line][itr2]
            if diff == True:
                ax = sns.lineplot(
                    x=range(len(numpy.diff(replicates_mean[line]))),
                    y=numpy.diff(replicates_mean[line]),
                    label=line,
                )
            else:
                print(line)

                line_l = line.replace("line_margin_", "Margin Length: ")
                """
                sns.lineplot(
                    x=range(len(replicates_mean[line])),
                    y=replicates_mean[line],
                    label=line,
                )
                """
        plt.ylabel("Odds Ratio")
        plt.xlabel("Generation #")

        # plt.legend(loc="upper left")
        plt.title(file_n)
        plt.savefig(f"{file_n}.png", dpi=300)
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
