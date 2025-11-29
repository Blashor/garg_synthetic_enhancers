import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy

super_ce = {"x": [], "y": []}
stitch_ce = {"x": [], "y": []}


def main():
    ce_obj = []
    clusters = []
    cluster_split = 200
    starting_stat = -1
    custom_labels = []
    with open("../crispr/super_info/ce_rank.txt") as file:
        for line in file:
            cell = line.strip().split("\t")
            start_ends = []
            if cell[-1] == "Super":
                for ce in cell[0:2]:
                    start, end = map(int, ce.split(":")[1].split("-"))
                    ce_len = end - start
                    start_ends.append((start, end))
                    # ce_score = float(cell[3])
                    # ce_obj.append((ce_score, ce_len))
                start_ends = sorted(start_ends, key=lambda item: item[0], reverse=False)
                print(start_ends)
                dist = start_ends[1][0] - start_ends[0][1]
                ce_obj.append((float(cell[3]), dist))

    ce_obj = sorted(ce_obj, key=lambda item: item[1], reverse=True)
    for ce_score, ce_len in ce_obj:
        if starting_stat == -1:
            starting_stat = ce_len
            custom_labels.append(str(starting_stat))
            cluster = []
        # if starting_stat - ce_len > cluster_split:
        if len(cluster) == 150:
            starting_stat = ce_len
            print(starting_stat)
            custom_labels.append(str(starting_stat))
            clusters.append(cluster)
            cluster = []
        else:
            cluster.append(ce_score)
    clusters.append(cluster)
    print("Plot")

    for cluster in clusters:
        print(len(cluster), numpy.mean(cluster))
        # custom_labels.append(str(round(numpy.mean(cluster))))

    sns.boxplot(data=clusters, width=0.5, linewidth=0.5)
    plt.xticks(range(len(custom_labels)), custom_labels)
    plt.show()


main()
