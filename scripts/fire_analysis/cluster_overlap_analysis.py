import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#
# Generate table for osn and fire overlaps
#


def chia_to_bed():
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


itr = 0
OSN_file = "/Users/blake/Documents/gargLab/FIBERseq/fiberSeqBeds/Analysis/youngSEportedmm10.bed"
TAD_file = "data/dowen_chia_mm10lift.bed"
cluster_file = "../metagene/fire_ce.svg3.bed"
# cluster_file = "/Users/blake/Downloads/fire_ry_comparable.svg3.bed"

os.system(f"bedtools intersect -a {OSN_file} -b {cluster_file} -wao > data/cluster_osn_overlap.bed")
os.system(f"bedtools intersect -a {TAD_file} -b {cluster_file} -wao > data/cluster_tad_overlap.bed")


# cluster_overlap_analysis
def osn_bar():
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
        label = ",".join(inters)
        if label == ".":
            label = "N/A"
        if label not in count_osn:
            count_osn[label] = 0
        count_osn[label] += 1

    count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
    keys = list(count_osn.keys())
    values = list(count_osn.values())
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


osn_supers = {}
with open("data/cluster_tad_overlap.bed") as file:
    for line in file:
        print(line)
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
    label = "_".join(inters)
    if label not in count_osn:
        count_osn[label] = 0
    count_osn[label] += 1

count_osn = dict(sorted(count_osn.items(), key=lambda item: item[1], reverse=True))
keys = list(count_osn.keys())
values = list(count_osn.values())
# /Users/blake/Documents/gargLab/FIRE_story_v1/data/cluster_rank.txt
# Plot bar chart
plt.bar(keys, values)

# Add numbers on top of the bars
for i in range(len(keys)):
    plt.text(keys[i], values[i], str(values[i]), ha="center", va="bottom")

# Add labels and title
plt.xlabel("Fire Clusters")
plt.title("CHIA Intersects")

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show plot
plt.tight_layout()
plt.show()

osn_bar()
