from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import os
import seaborn as sns
import scipy
import math


def ce_rank_to_bed(file_name, ignore_super=False):
    super_set = set()
    super_set_filter = set()
    if ignore_super == False:
        super_filter = "Super"
    else:
        super_filter = ""
    with open(file_name) as file:
        for line in file:
            cell = line.split("\t")
            if super_filter in line:
                # print(cell, cell[2])
                chrom, c_range = cell[2].split(":")
                start, end = map(int, c_range.split("-"))
                if f"{chrom}\t{start}\t{end}" not in super_set_filter:
                    if float(cell[-3]) != 0:
                        # super_set.add(f"{chrom}\t{start}\t{end}\t{cell[-3]}\n")
                        super_set.add(f"{chrom}\t{start}\t{end}\t{math.log(float(cell[-3]),2)}\n")
                    else:
                        pass
                        # super_set.add(f"{chrom}\t{start}\t{end}\t{0}\n")
                    super_set_filter.add(f"{chrom}\t{start}\t{end}")
                else:
                    pass
                    # print("000ooo000")
    print(len(list(super_set_filter)))
    with open(f"{file_name}.bed", "w") as file:
        file.writelines(list(super_set))
    return len(super_set)


def get_file_len(file_name):
    i = 0
    with open(file_name) as file:
        for line in file:
            if "\t" in line:
                i += 1
    return i


def main(f1="", f2="", output="", ignore_super=False):
    file1 = f1
    file2 = f2
    plt.figure()
    set1 = set()
    set2 = set()

    file1_len = ce_rank_to_bed(file1, ignore_super)

    file2_len = ce_rank_to_bed(file2, ignore_super)
    print(file1_len, file2_len)
    os.system(f"bedtools intersect -a {file1}.bed -b {file2}.bed -wao > ../data/intersect.bed")
    intersect_count = 0
    x = []
    y = []

    with open("../data/intersect.bed") as file:
        for line in file:
            print(line)
            cell = line.split("\t")
            if cell[3] == ".":
                cell[3] = 0
            if cell[7] == ".":
                cell[7] = 0
            x.append(float(cell[3]))
            y.append(float(cell[7]))
            print(x[-1], y[-1])
    sns.scatterplot(
        x=x,
        y=y,
        s=2,
        alpha=0.25,
        color=(157 / 255, 36 / 255, 36 / 255),
        edgecolor=(157 / 255, 36 / 255, 36 / 255),
    )

    ax = sns.regplot(
        x=x,
        y=y,
        color=None,
        order=1,
        scatter_kws={"color": None, "s": 10, "alpha": 0},
        line_kws={"linestyle": "--", "color": (157 / 255, 36 / 255, 36 / 255)},
    )
    print(len(x))
    print(scipy.stats.linregress(x=x, y=y, alternative="two-sided"))
    plt.show()


"/Users/blake/Documents/gargLab/Figures_for_Box/Figure_1_Ranking/data/"
"""
main_venn3(
    f1="../data/3runs-FDR-FIRE-peaks.bed",
    f2="../data/SRX23682289_atac.05.bed",
    f3="../data/SRX10040677_danse.05.bed",
    output="venn_fire_atac.svg",
)
"""
"""
main(
    f1="../data/GM12878_ce_rank.txt",
    f2="../data/K562_ce_rank.txt",
    output="../figures/s5c_hg38_not_super.svg",
    ignore_super=False,
)
"""
# main(f1="../data/3runs_ce_rank.txt", f2="../data/ce_rank_no_dist_correctionyoungSE_OSN.txt", ignore_super=False)
main(
    f1="../data/3runs_ce_rank.txt",
    f2="/Users/blake/Documents/gargLab/for_github/ranking_fisher_1/ranked_pipeline/data/ce_rankrun4_jun27.txt",
    # f2="../data/ce_rank_michele_compatible.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=True,
)
