import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy as np
import Bio.SeqUtils
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import seaborn as sns


def smoothed_gradient(y, x=None, span=20):
    """
    Compute gradient using points `span` indices ahead and behind.
    Essentially (f(x+h) - f(x-h)) / (2h)

    Parameters:
        y: array-like, data values
        x: array-like or None, optional x-values
        span: int, number of points away on either side to use for gradient

    Returns:
        grad: numpy array, smoothed gradient
    """
    y = np.asarray(y)
    if x is None:
        x = np.arange(len(y))
    else:
        x = np.asarray(x)

    grad = np.empty_like(y)
    grad[:] = np.nan  # edges will be undefined

    for i in range(span, len(y) - span):
        dy = y[i + span] - y[i - span]
        dx = x[i + span] - x[i - span]
        grad[i] = dy / dx

    return grad


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
    rep_holder = {}
    for gen_num, gen in enumerate(gen_holder["score"][:300]):
        for rep_num, reps in enumerate(gen):
            # print(gen_num, rep_num)
            if rep_num not in rep_holder:
                rep_holder[rep_num] = []
            rep_holder[rep_num].append(np.sqrt(np.average(reps)))
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 8), sharex=True)
    for ax in axes:
        # Hide top and right spines
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        # Style bottom and left spines
        ax.spines["bottom"].set_linewidth(1.5)
        ax.spines["left"].set_linewidth(1.5)

        # Style tick marks
        ax.tick_params(axis="both", width=2)
    for rep in rep_holder:
        x = range(len(rep_holder[rep]))
        y = rep_holder[rep]
        # dy = smoothed_gradient(rep_holder[rep], x)
        # d2y = smoothed_gradient(dy, x)
        smooth_y = gaussian_filter1d(y, sigma=10)
        dy = np.gradient(smooth_y, x)
        d2y = np.gradient(dy, x)
        """
        sns.lineplot(
            x=x,
            y=d2y,
            alpha=0.8
            # label=line,
        )
        """
        sns.lineplot(x=x, y=smooth_y, ax=axes[0], alpha=0.8)
        axes[0].set_ylabel("Odds Ratio")

        # Plot dy
        sns.lineplot(x=x, y=dy, ax=axes[1], alpha=0.8)
        axes[1].set_ylabel("First Derivative")

        # Plot d2y
        sns.lineplot(x=x[:298], y=d2y[:298], ax=axes[2], alpha=0.8)
        axes[2].set_ylabel("Second Derivative")
        axes[2].set_xlabel("Generation")

        plt.tight_layout()

    plt.show()
    plt.savefig(output)


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
        sns.lineplot(
            x=range(len(rep_holder[rep])),
            y=rep_holder[rep],
            alpha=0.8
            # label=line,
        )
    plt.savefig("../figures/3e_gc_per_gen.svg")


def main(file_n):
    gen_add = 0
    surv = 4
    # file_n = "INSIDE_jun13.txt"
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
                gen, line_rep, seq, score = line.split("\t")
                gen = int(gen.split(" ")[-1])
                if "a" in (line_rep[-1]) or "b" in (line_rep[-1]) or "c" in (line_rep[-1]):
                    line_rep = str(line_rep[:-1])
                if line_rep.split("rep")[-1] not in rep_tracker:
                    # print(line_rep)
                    rep_tracker[line_rep.split("rep")[-1]] = rep_number
                    rep_number += 1
                rep = rep_tracker[line_rep.split("rep")[-1]]
                score = float(score)
                if last_gen != gen:
                    rep_holder1 = [[] for x in range(40)]
                    rep_holder2 = [[] for x in range(40)]
                    last_gen = gen
                    real_gen += 1
                    gen_holder["score"].append(rep_holder1)
                    gen_holder["seq"].append(rep_holder2)
                    # print(real_gen)
                print(rep)
                gen_holder["seq"][real_gen][rep - 1].append(seq)
                gen_holder["score"][real_gen][rep - 1].append(score)
    # print(gen_holder["score"])
    return gen_holder


gen_holder = main("../data/INSIDE_jun17.txt")

# calc_gc(gen_holder)
score_by_gen(gen_holder, "../figures/s7c_inflection.svg")
