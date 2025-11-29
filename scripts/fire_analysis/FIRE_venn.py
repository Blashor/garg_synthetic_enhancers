from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import os


def ce_rank_to_bed(file_name, ignore_super=False):
    super_set = set()
    super_set_filter = set()
    with open(file_name) as file:
        for line in file:
            cell = line.split("\t")

            # print(cell, cell[2])
            chrom, c_range = cell[2].split(":")
            start, end = map(int, c_range.split("-"))
            # print(cell[-1])
            if f"{chrom}\t{start}\t{end}" not in super_set_filter:
                if cell[-1] != "Super\n":
                    cell[-1] = cell[-3]
                super_set.add(f"{chrom}\t{start}\t{end}\t{cell[-1]}\n")
                super_set_filter.add(f"{chrom}\t{start}\t{end}")
            else:
                pass

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

    file1_len = ce_rank_to_bed(file1)

    file2_len = ce_rank_to_bed(file2, ignore_super)
    print(file1_len, file2_len)
    os.system(f"bedtools intersect -a {file1}.bed -b {file2}.bed > ../data/intersect.bed")
    intersect_count = 0
    with open("../data/intersect.bed") as file:
        for line in file:
            if "chr" in line:
                intersect_count += 1
    for i in range(intersect_count):
        set1.add(i)
        set2.add(i)
    for i in range(file1_len - intersect_count):
        set1.add(i + 1000000)
    for i in range(file2_len - intersect_count):
        set2.add(i + 5000000)
    label1 = file1.split("/")[-1].split(".txt")[0]
    label2 = file2.split("/")[-1].split(".txt")[0]
    venn2([set1, set2], set_labels=(label1, label2))
    plt.savefig(output)
    plt.show()


def main2(f1="", f2="", output="", ignore_super=False):
    file1 = f1
    file2 = f2
    plt.figure()
    set1 = set()
    set2 = set()

    file1_len = ce_rank_to_bed(file1)

    file2_len = ce_rank_to_bed(file2, ignore_super)
    # print(file1_len, file2_len)
    os.system(f"bedtools intersect -a {file1}.bed -b {file2}.bed -wao > ../data/intersect.bed")
    intersect_count = 0
    file1_len = 0
    file2_len = 0
    with open("../data/intersect.bed") as file:
        for line in file:
            cell = line.split("\t")
            if cell[3] == "Super":
                if cell[-2] == "Super":
                    intersect_count += 1
                    # print("both_super")
                else:
                    file1_len += 1
            elif cell[-2] == "Super":
                file2_len += 1

    print(intersect_count)
    for i in range(intersect_count):
        set1.add(i)
        set2.add(i)
    for i in range(file1_len - intersect_count):
        set1.add(i + 1000000)
    for i in range(file2_len - intersect_count):
        set2.add(i + 5000000)
    label1 = file1.split("/")[-1].split(".txt")[0]
    label2 = file2.split("/")[-1].split(".txt")[0]
    venn2([set1, set2], set_labels=(label1, label2))
    plt.savefig(output)
    plt.show()


def main_venn3(f1="", f2="", f3="", output="", ignore_super=False):
    file1 = f2
    file2 = f3
    plt.figure()
    set1 = set()
    set2 = set()

    file1_len = get_file_len(file1)

    file2_len = get_file_len(file2)
    print(file1_len, file2_len)
    for file1, file2 in ((file1, file2), (file2, file1)):
        os.system(f"bedtools intersect -a {file1} -b {file2} -u > ../data/intersect_u.bed")
        os.system(f"bedtools intersect -a {file1} -b {file2} -v > ../data/intersect_v.bed")
        print("u", get_file_len("../data/intersect_u.bed"), "v", get_file_len("../data/intersect_v.bed"))
        print()
    intersect_count = 0
    with open("../data/intersect.bed") as file:
        for line in file:
            if "chr" in line:
                intersect_count += 1
    for i in range(intersect_count):
        set1.add(i)
        set2.add(i)
    for i in range(file1_len - intersect_count):
        set1.add(i + 1000000)
    for i in range(file2_len - intersect_count):
        set2.add(i + 5000000)
    label1 = file1.split("/")[-1].split(".txt")[0]
    label2 = file2.split("/")[-1].split(".txt")[0]
    venn2([set1, set2], set_labels=(label1, label2))
    plt.savefig(output)
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
main2(
    f1="../data/3runs_ce_rank.txt",
    f2="/Users/blake/Downloads/ce_rank_sep1.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=False,
)

main2(
    f1="../data/ce_rank_michele_compatible.txt",
    f2="/Users/blake/Downloads/ce_rank_sep1.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=False,
)
main2(
    f1="../data/3runs_ce_rank.txt",
    f2="../data/ce_rank_michele_compatible.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=False,
)
main2(f1="../data/GM12878_ce_rank.txt", f2="../data/K562_ce_rank.txt", output="../figures/s5c_hg38_not_super.svg")
# main(f1="../data/3runs_ce_rank.txt", f2="../data/ce_rank_no_dist_correctionyoungSE_OSN.txt", ignore_super=False)
main2(
    f1="../data/3runs_ce_rank.txt",
    f2="../data/ce_rank_michele_compatible.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=False,
)


main2(
    f1="../data/ce_rank_michele_compatible.txt",
    f2="/Users/blake/Documents/gargLab/for_github/ranking_fisher_1/ranked_pipeline/data/ce_rankrun4_jun27.txt",
    output="../figures/s5b_mm10.svg",
    ignore_super=True,
)
