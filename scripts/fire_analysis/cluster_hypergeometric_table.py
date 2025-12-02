import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# hypergeometric


def FIRE_super_to_bed():
    bed_lines = set()
    with open("/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/3runs_ce_rank.txt") as file:
        for line in file:
            if "Super" in line:
                cell = line.strip().split("\t")
                chrom, c_range = cell[0].split(":")
                start, end = c_range.split("-")
                bed_lines.add(f"{chrom}\t{start}\t{end}\n")
                #
                chrom, c_range = cell[1].split(":")
                start, end = c_range.split("-")
                bed_lines.add(f"{chrom}\t{start}\t{end}\n")
    print(len(list(bed_lines)))
    with open("data/FIRE_super.bed", "w") as file:
        file.writelines(list(bed_lines))


def chia_to_bed():
    # Have to convert to mm10 after
    chia_lines = []
    with open("data/CHIA_pet_dowen_sup1.txt") as file:
        for line in file:
            cell = line.strip().split("\t")
            print("\t".join(cell[:3]))
            print("\t".join(cell[3:]))
            chia_lines.append("\t".join(cell[:3]) + "\n")
            chia_lines.append("\t".join(cell[3:]) + "\n")
    with open("data/CHIA.bed", "w") as file:
        file.writelines(chia_lines)


def osn_hyp():
    osn_supers = {}
    with open("data/cluster_osn_overlap.bed") as file:
        for line in file:
            cell = line.split("\t")
            osn = "\t".join(cell[:3])
            cluster = cell[-2].split("_")[-1]
            if osn not in osn_supers:
                osn_supers[osn] = set([])
            osn_supers[osn].add(cluster)
    count_osn = {}
    for osn in osn_supers:
        inters = list(osn_supers[osn])
        inters = sorted(inters)
        for label in inters:
            if label not in count_osn:
                count_osn[label] = 0
            count_osn[label] += 1

    count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
    print("RY", count_osn)


def osn_hyp():
    osn_supers = {}
    with open("data/cluster_osn_overlap.bed") as file:
        for line in file:
            cell = line.split("\t")
            osn = "\t".join(cell[:3])
            cluster = cell[-2].split("_")[-1]

            if osn not in osn_supers:
                osn_supers[osn] = set([])
            osn_supers[osn].add(cluster)
    count_osn = {}
    cluster_overlap_count = {}
    for osn in osn_supers:
        inters = list(osn_supers[osn])
        inters = sorted(inters)
        cluster_name = "_".join(inters)
        if cluster_name == ".":
            cluster_name = "N/A"
        if cluster_name not in cluster_overlap_count:
            cluster_overlap_count[cluster_name] = 0
        cluster_overlap_count[cluster_name] += 1
        for label in inters:
            if label not in count_osn:
                count_osn[label] = 0
            count_osn[label] += 1

    count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
    print("RY", count_osn)

    cluster_overlap_count = dict(sorted(cluster_overlap_count.items(), key=lambda item: item[1], reverse=True))
    keys = list(cluster_overlap_count.keys())
    values = list(cluster_overlap_count.values())
    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_subplot(111)  # 111 means 1 row, 1 column, 1st subplot
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)  # Bottom (x-axis) spine
    ax.spines["left"].set_linewidth(1.5)  # Left (y-axis) spine
    ax.tick_params(axis="both", width=2)

    # Plot bar chart
    plt.bar(keys, values, color="black")

    # Add numbers on top of the bars
    for i in range(len(keys)):
        plt.text(keys[i], values[i], str(values[i]), ha="center", va="bottom")

    # Add labels and title
    plt.xlabel("Fire Clusters")
    plt.title("OSN Super Intersects")

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=-45)

    # Show plot
    plt.tight_layout()
    plt.show()


# osn_bar()


def chia_hyp():
    osn_supers = {}
    with open("data/cluster_tad_overlap.bed") as file:
        for line in file:
            cell = line.split("\t")
            osn = "\t".join(cell[:3])
            cluster = cell[-2].split("_")[-1]
            # print(cell)
            if osn not in osn_supers:
                osn_supers[osn] = set([])
            osn_supers[osn].add(cluster)
    count_osn = {}
    for osn in osn_supers:
        inters = list(osn_supers[osn])
        inters = sorted(inters)
        for label in inters:
            if label not in count_osn:
                count_osn[label] = 0
            count_osn[label] += 1

    count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
    print("TAD", count_osn)


def FIRE_SUPER_hyp():
    osn_supers = {}
    with open("data/cluster_FIRE_SE_overlap.bed") as file:
        for line in file:
            cell = line.split("\t")
            osn = "\t".join(cell[:3])
            cluster = cell[-2].split("_")[-1]
            if osn not in osn_supers:
                osn_supers[osn] = set([])
            osn_supers[osn].add(cluster)
    count_osn = {}
    for osn in osn_supers:
        inters = list(osn_supers[osn])
        inters = sorted(inters)
        for label in inters:
            if label not in count_osn:
                count_osn[label] = 0
            count_osn[label] += 1

    count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
    print("FIRE SE", count_osn)


FIRE_super_to_bed()
chia_to_bed()

itr = 0
OSN_file = "/Users/blake/Documents/gargLab/FIBERseq/fiberSeqBeds/Analysis/youngSEportedmm10.bed"
FIRE_SE_file = "data/FIRE_super.bed"
TAD_file = "data/dowen_chia_mm10lift.bed"


# OSN_file = "/Users/blake/Documents/gargLab/FIBERseq/fiberSeqBeds/Analysis/youngSEportedmm10.bed"
# TAD_file = "data/dowen_chia_mm10lift.bed"
# cluster_file = "../metagene/fire_ce.svg3.bed"  # red these fires use our filters for downstream analysis
# cluster_file = "/Users/blake/Downloads/fire_ry_comparable.svg3.bed"


def filter_deeptools_cluster_file(cluster_file):
    filter_set = set()
    fire_intergenic = "/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/3runV2_intergenic.bed"
    with open(fire_intergenic) as file:
        for line in file:
            # print(line)
            cell = line.strip().split("\t")
            filter_set.add(f"{cell[0]}:{cell[1]}-{cell[2]}")
    lines_to_add = []
    with open(cluster_file) as file:
        for line in file:
            line = line.strip()
            cell = line.split("\t")
            if cell[3] in filter_set:
                lines_to_add.append(line + "\n")
            else:
                print(line)
    print(len(lines_to_add))
    with open("filter_cluster.txt", "w") as file:
        file.writelines(lines_to_add)


pre_cluster_file = "/Users/blake/Documents/gargLab/fire_story_v3/fire_3runs_BOX.svg3.bed"
filter_deeptools_cluster_file(pre_cluster_file)
cluster_file = "filter_cluster.txt"

clusters_count = {}
with open(cluster_file) as file:
    for line in file:
        cell = line.split("\t")
        if cell[-1] not in clusters_count:
            clusters_count[cell[-1]] = 0
        clusters_count[cell[-1]] += 1


os.system(f"bedtools intersect -a {OSN_file} -b {cluster_file} -wao > data/cluster_osn_overlap.bed")
os.system(f"bedtools intersect -a {FIRE_SE_file} -b {cluster_file} -wao > data/cluster_FIRE_SE_overlap.bed")
os.system(f"bedtools intersect -a {TAD_file} -b {cluster_file} -wao > data/cluster_tad_overlap.bed")

itr = 0
with open(FIRE_SE_file) as file:
    for line in file:
        itr += 1
print(itr)
itr = 0
with open(TAD_file) as file:
    for line in file:
        itr += 1
print(itr)
itr = 0
with open(OSN_file) as file:
    for line in file:
        itr += 1
print(itr)
clus = {}
with open(cluster_file) as file:
    for line in file:
        c = line.strip()[-1]

        if c not in clus:
            clus[c] = 0
        clus[c] += 1
print(clus)


FIRE_SUPER_hyp()
osn_hyp()
chia_hyp()
print(clusters_count)
