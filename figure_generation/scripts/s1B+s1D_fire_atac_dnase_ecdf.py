import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import mpl_style
import numpy as np


def main(file_name="../data/3runs-FIRE.bed", file_out=""):
    mpl_style.set_style(figsize=(4.5, 4.5))
    # plt.figure()
    """
    size_arr = []
    with open("../data/fire-fibers.bed") as file:
        for line in file:
            cell = line.split("\t")
            size_arr.append(int(cell[2]) - int(cell[1]))
    sns.ecdfplot(data=size_arr)
    plt.show()
    """
    size_arr = []
    with open(file_name) as file:
        for line in file:
            cell = line.split("\t")
            size_arr.append(int(cell[2]) - int(cell[1]))
    print(np.median(size_arr), np.quantile(size_arr, [0.10, 0.50, 0.90]))
    sns.ecdfplot(data=size_arr, color="black")
    plt.ylim([0, 1.1])
    plt.xlim([0, 30000])
    plt.savefig(file_out)
    plt.show()


def mains(
    file_names=[
        "../data/3runs-FDR-FIRE-peaks.bed",
        "../data/SRX10040677.05.bed",
        "../data/SRX23682289.05.bed",
    ],
    file_out="",
):
    mpl_style.set_style(figsize=(4.5, 4.5))
    # plt.figure()
    """
    size_arr = []
    with open("../data/fire-fibers.bed") as file:
        for line in file:
            cell = line.split("\t")
            size_arr.append(int(cell[2]) - int(cell[1]))
    sns.ecdfplot(data=size_arr)
    plt.show()
    """

    for file_name in file_names:
        size_arr = []
        with open(file_name) as file:
            for line in file:
                cell = line.split("\t")
                try:
                    size_arr.append(int(cell[2]) - int(cell[1]))
                except Exception as e:
                    pass
                # size_arr.append(int(cell[2]) - int(cell[1]))
        print(np.median(size_arr), np.quantile(size_arr, [0.10, 0.50, 0.90]))
        if "FIRE" in file_name:
            label = "Fire"
            sns.ecdfplot(data=size_arr, color="#9d1e1eff", label=f"{label} Median: {int(np.median(size_arr))} bp")
        if "dnase" in file_name:
            label = "DNase"
        if "atac" in file_name:
            label = "ATAC"
        if "FIRE" not in file_name:
            sns.ecdfplot(data=size_arr, label=f"{label} Median: {int(np.median(size_arr))} bp")
    plt.ylim([0, 1.1])
    plt.xlim([0, 1500])
    plt.legend()
    plt.show()
    plt.savefig(file_out)


mains(file_out="../figures/s2d_fire_atac_dnase.svg")
# sys.exit()
main(file_name="../data/3runs-fire-fibers.bed", file_out="../figures/s2d_fire_fibers_ecdf.svg")
# main(file_out="../figures/s2d_fire_region.svg")
