import math
import sys
import json
import scipy
import numpy
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import fastcluster
from sklearn.cluster import KMeans
import pandas as pd
from scipy.stats.stats import pearsonr
import random
from scipy.optimize import curve_fit


#
# Plot Hockey Stick from ranking_pipeline ranked ce pair .txt output
#


def numPts_below_line(myVector, slope, x):
    yPt = myVector[x]
    b = yPt - (slope * x)
    xPts = numpy.arange(1, len(myVector) + 1)
    return numpy.sum(myVector <= (xPts * slope + b))


def calculate_cutoff(inputVector, drawPlot=True, return_y_thresh=False, savePlot="", drawTangent=True):
    inputVector = numpy.sort(inputVector)
    inputVector[inputVector < 0] = 0

    slope = (numpy.max(inputVector) - numpy.min(inputVector)) / len(inputVector)

    result = numpy.vectorize(lambda x: numPts_below_line(inputVector, slope, x))(numpy.arange(len(inputVector)))
    xPt = numpy.floor(numpy.argmin(result)) + 1
    y_cutoff = inputVector[int(xPt - 1)]

    if drawPlot:
        # figsize=(3, 3)
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
        ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
        ax.tick_params(axis="both", width=2)
        if drawTangent:
            b = y_cutoff - (slope * xPt)
            x_vals = numpy.array([0, xPt + 10000])  # X range for the slope line
            y_vals = slope * x_vals + (b)  # Corresponding Y values
            ax.plot(x_vals, y_vals, "k--", linewidth=2, label="Slope Line", alpha=0.5)
            # Plot inputVector
        ax.plot(numpy.arange(1, len(inputVector) + 1), inputVector, "k-", linewidth=2, solid_capstyle="round", label="")

        # Plot supers
        supers = inputVector[inputVector > y_cutoff]
        ax.plot(
            numpy.arange(xPt, len(supers) + xPt),
            supers,
            "r-",
            linewidth=2.1,
            solid_capstyle="round",
            label="Super Co-accessible",
        )
        # plt.xlim([-100, len(inputVector) + 1000])
        # plt.ylim([-100, inputVector[-1]])

        plt.xlabel("Co-Accessible Regions Ranked By Odds Ratio")
        plt.ylabel("Odds Ratio")
        plt.legend()

        if savePlot:
            plt.tight_layout()
            plt.savefig(savePlot)
            plt.show()

    if return_y_thresh:
        return y_cutoff
    return xPt


def main(fp="../data/cluster_rank.txt", save="../figures/clusters_ranked.svg", tangent=False):
    odds_ratios = []
    with open(fp) as file:
        for line in file:
            cell = line.split("\t")
            odds_ratios.append(float(cell[3]))
    plt.figure()
    y_cutoff = calculate_cutoff(
        odds_ratios,
        drawPlot=True,
        return_y_thresh=True,
        savePlot=save,
        drawTangent=tangent,
    )


# main(fp="../data/3runs_ce_rank.txt", save="../figures/1c_3runs_ce_rank.svg")
# main(fp="../data/GM12878_ce_rank.txt", save="../figures/1e_GM12878_ce_rank.svg")

# 1f: K562 Hockey green
# main(fp="../data/K562_ce_rank.txt", save="../figures/1f_K562_ce_rank.svg")
# print(sns.cm.rocket)
# main()
