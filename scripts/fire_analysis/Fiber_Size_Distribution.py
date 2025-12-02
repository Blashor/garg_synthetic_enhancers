import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import mpl_style
import numpy as np


#
# eCDF of bed file element sizes
#
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
    plt.xlim([0, 1500])
    plt.savefig(file_out)


def mains(
    file_names=[
        "../data/3runs-FDR-FIRE-peaks.bed",
        "/Users/blake/Downloads/SRX10040677_dnase.05.bed",
        "/Users/blake/Downloads/SRX23682289_atac2.05.bed",
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
sys.exit()
# main(file_name="../data/3runs-fire-fibers.bed", file_out="../figures/s2d_fire_fibers_ecdf.svg")
main(file_out="../figures/s2d_fire_region.svg")

# sys.exit()
genome_size = {}
with open("/Users/blake/Downloads/GRCm38.primary_assembly.genome.fa.fai") as file:
    for line in file:
        cell = line.split("\t")
        genome_size[cell[0]] = int(cell[1])
total = 0
""
totals = {}
median_count = {}
with open("../data/fiber_cov.bg") as file:
    for line in file:
        if "chr" in line:
            cell = line.strip().split("\t")
            if len(cell) == 4:
                if cell[0] not in totals:
                    totals[cell[0]] = 0
                #
                length = int(cell[2]) - int(cell[1])
                totals[cell[0]] += length
                #
                cov = int(cell[-1])
                if cov not in median_count:
                    median_count[cov] = 0
                median_count[cov] += length
                # total += cov
# fill in zeros
for chrom in totals:
    if 0 not in median_count:
        median_count[0] = 0
    else:
        median_count[0] += genome_size[chrom] - totals[chrom]

sorted_counts = sorted(median_count.items())  # Sort coverage levels
total_length = sum(median_count.values())  # Total number of base pairs
half_length = total_length / 2

cumulative = 0
median = None

for cov, length in sorted_counts:
    cumulative += length
    if cumulative >= half_length:
        median = cov
        break
print(sorted_counts)
print("Median Coverage:", median)
print(total)
print(total / 2730871774)
#
#
