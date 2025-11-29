import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import json
import sys


def main():
    ce_stitch_sizes = []
    ce_stitch_dist = []
    ce_super_sizes = []
    ce_super_dist = []
    with open("../FIRE_story_v1/data/ce_rank.txt") as file:
        for line in file:
            cell = line.split("\t")
            start1, end1 = map(int, cell[0].split(":")[1].split("-"))
            start2, end2 = map(int, cell[1].split(":")[1].split("-"))
            fire1_size = end1 - start1
            fire2_size = end2 - start2
            dist_between = start2 - end1
            if dist_between > 0:
                if cell[5].strip() == "Super":
                    ce_super_sizes.append(fire1_size)
                    ce_super_sizes.append(fire2_size)
                    ce_super_dist.append(dist_between)
                else:
                    ce_stitch_sizes.append(fire1_size)
                    ce_stitch_sizes.append(fire2_size)
                    ce_stitch_dist.append(dist_between)

    mean = np.mean(ce_stitch_sizes)
    q25, q50, q75 = np.percentile(ce_stitch_sizes, [25, 50, 75])
    print(q25, q50, q75, mean)

    mean = np.mean(ce_stitch_dist)
    q25, q50, q75 = np.percentile(ce_stitch_dist, [25, 50, 75])
    print(q25, q50, q75, mean)
    print("Supers:")

    mean = np.mean(ce_super_sizes)
    q25, q50, q75 = np.percentile(ce_super_sizes, [25, 50, 75])
    print(q25, q50, q75, mean)

    mean = np.mean(ce_super_dist)
    q25, q50, q75 = np.percentile(ce_super_dist, [25, 50, 75])
    print(q25, q50, q75, mean)
    sns.ecdfplot(ce_stitch_dist)
    sns.ecdfplot(ce_super_dist)
    plt.show()


main()
