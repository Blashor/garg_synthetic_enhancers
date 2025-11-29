import math
import sys
import json
import scipy
import numpy
from matplotlib import pyplot as plt
import seaborn as sns
import mpl_style
import numpy as np


def jitter_array(arr, jitter_strength=0.03):
    arr = np.array(arr)
    log_arr = np.log10(arr + 1e-8)  # avoid log(0)
    jitter = np.random.normal(loc=0, scale=jitter_strength, size=arr.shape)
    log_arr_jittered = log_arr + jitter
    jittered_arr = 10**log_arr_jittered
    return np.clip(jittered_arr, a_min=1e-3, a_max=1e3)


def main(
    fps=["../data/3runs_ce_rank.txt"], labels=[], colors=False, save="../figures/clusters_ranked.svg", supers_only=False
):
    plt.figure()
    mpl_style.set_style((4, 4))
    if supers_only == False:
        supers_only = [False] * len(fps)
    if colors == False:
        colors = [False] * len(fps)
    for fp, label, c, s in zip(fps, labels, colors, supers_only):
        odds_ratios = []
        with open(fp) as file:
            for line in file:
                cell = line.split("\t")
                if float(cell[3]) != 0:
                    odds_ratios.append(float(cell[3]))
        odds_ratios = jitter_array(odds_ratios)
        if c != False:
            sns.ecdfplot(data=odds_ratios, label=label, color=c)
        else:
            sns.ecdfplot(data=odds_ratios, label=label)
        plt.xlabel("Odds Ratio")
        plt.ylim([0, 1.1])
        plt.xlim([10**-2, 10**2])
        plt.xscale("log")
    plt.legend()
    plt.show()
    plt.savefig(save)


fps = [
    "../data/3runs_ce_rank_no_dist_correction.txt",
    "/Users/blake/Documents/gargLab/for_github/ranking_fisher_1/ranked_pipeline/data/ce_rank_no_dist_correctionfire_stitched_control0.txt",
]
main(fps, labels=["Stitched", "Control"], save="../figures/s3d_fire_ecdf.svg")
