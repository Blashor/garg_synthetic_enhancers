import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def set_style(figsize=(7, 6)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)


itr = 0
OSN_file = "/Users/blake/Documents/gargLab/FIBERseq/fiberSeqBeds/Analysis/youngSEportedmm10.bed"
cluster_file = "/Users/blake/Downloads/fire_ry_comparable.svg3.bed"


os.system(f"bedtools intersect -a {OSN_file} -b {cluster_file} -wao > data/cluster_osn_overlap.bed")


def osn_bar():
    i = 0
    RY_supers_to_FIRE = {}
    osn_supers = {}
    with open("data/cluster_osn_overlap.bed") as file:
        for line in file:
            cell = line.split("\t")
            osn = "\t".join(cell[:3])
            cluster = cell[-2].split("_")[-1]
            i += 1
            # print(cell)
            # print(i, cell[6], cluster)
            RY_supers_to_FIRE[cell[6]] = (cluster, osn)
    x_ry = []
    y_ry = []
    y_all = []
    y_super = []
    s_i = 0
    osn_used = set([])
    with open("/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/3runs_ce_rank.txt") as file:
        for line in file:
            if ":" in line:
                ceA, ceB, stitch, OR, rank, isSuper = line.split("\t")
                OR = float(OR)
                y_all.append(float(OR))

                if ceA in RY_supers_to_FIRE:
                    osn = RY_supers_to_FIRE[ceA][1]
                    if osn not in osn_used:
                        osn_used.add(osn)
                        x_ry.append(int(rank))
                        y_ry.append(float(OR))
                        s_i += 1
                        print(s_i, int(rank) / 100000, RY_supers_to_FIRE[ceA][0], OR)
                        del RY_supers_to_FIRE[ceA]
                if ceB in RY_supers_to_FIRE:
                    osn = RY_supers_to_FIRE[ceB][1]
                    if osn not in osn_used:
                        osn_used.add(osn)
                        s_i += 1
                        print(s_i, int(rank) / 100000, RY_supers_to_FIRE[ceB][0], OR)
                        del RY_supers_to_FIRE[ceB]
                if isSuper.strip() != "Super":
                    isSuper = "w"
                else:
                    y_super.append(OR)
    set_style()
    plt.scatter(x_ry, y_ry, s=5, alpha=0.2)
    plt.gca().invert_xaxis()
    plt.show()
    plt.figure()
    set_style()
    sns.ecdfplot(data=y_ry)
    sns.ecdfplot(data=y_all)
    sns.ecdfplot(data=y_super)
    plt.xscale("log")
    plt.show()
    # c_id = "_".join(map(str, sorted([cluster_labels[ceA], cluster_labels[ceB]]))) + f"{isSuper.strip()}"
    # c_id_no_super = "_".join(map(str, sorted([cluster_labels[ceA], cluster_labels[ceB]])))
    # cell = line.strip().split("\t")


osn_bar()
