import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy as np

# import Bio.SeqUtils

import matplotlib.pyplot as plt
import seaborn as sns


def set_style(figsize=(4, 3)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


labels = {
    "line1": "SNP 0.4%",
    "line2": "SNP 1%",
    "line3": "1-3 Duplications",
    "line4": "1-3 Duplications + SNP 1%",
    "line5": "SNP 5%",
    "line6": "SNP 20%",
}


def score_by_gen(gen_holder, output):
    set_style()
    for rep in gen_holder:
        print(rep)
        sns.lineplot(
            x=range(len(gen_holder[rep][:150])),
            y=np.sqrt(gen_holder[rep][:150] / 3),
            alpha=0.8,
            label=rep,
        )
    plt.savefig(output)
    plt.show()


def calc_gc(gen_holder):
    set_style()
    rep_holder = {}
    for gen_num, gen in enumerate(gen_holder["seq"][:300]):
        for rep_num, reps in enumerate(gen):
            # print(gen_num, rep_num)
            if rep_num not in rep_holder:
                rep_holder[rep_num] = []
            gc_values = list(map(lambda seq: 100 * Bio.SeqUtils.gc_fraction(seq), reps))
            rep_holder[rep_num].append(np.average(gc_values))
    for rep in rep_holder:
        print(rep)
        sns.lineplot(
            x=range(len(rep_holder[rep])),
            y=rep_holder[rep] / 3,
            alpha=0.8
            # label=line,
        )
    plt.savefig("../figures/3e_gc_per_gen.svg")


def main(file_n):
    gen_add = 1
    surv = 4
    # file_n = "../data/INSIDE_jun13.txt"
    # file_n = "INSIDE_graded_on4_3model.txt"
    gen_holder = {"score": [], "seq": []}
    last_gen = 0
    last_rep = 0
    real_gen = -1  # to deal with an earlier adding issue
    rep_number = 0
    rep_tracker = {}
    with open(file_n) as file:
        for line in file:
            if "Gen:" in line:
                cell = line.split("\t")
                name = cell[1].split("rep")[0]
                print(cell[1].split("rep")[-1])
                g = int(cell[0].split("Gen: ")[-1]) - 1
                print(g)
                if name not in rep_tracker:
                    rep_tracker[name] = np.zeros([300])
                if g < 300:
                    rep_tracker[name][g] += float(cell[-1])
    print(rep_tracker)
    return rep_tracker


examples = [
    "../data/INSIDE_margin.txt",
    "../data/INSIDE_len.txt",
    "../data/INSIDE_cs.txt",
    "../data/INSIDE_mut_mar26.txt",
]
"""
examples = [
    "../data/INSIDE_genomeA.out",
    "../data/INSIDE_genomeB.out",
    "../data/INSIDE_genomeC.out",
]
"""

for example in examples:
    print(example)
    gen_holder = main(example)
    name = example.split("/")[-1] + ".svg"
    # score_by_gen(gen_holder, f"../figures/{name}")
    score_by_gen(gen_holder, f"../figures/s6_{name}")
