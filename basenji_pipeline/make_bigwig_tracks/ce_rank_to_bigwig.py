# will start with bedgraph to see if that works with basenji
import os
import sys
import numpy


def main(file_path):
    ce_obj = {}
    just_path, just_file = os.path.split(file_path)
    print(just_path)
    bw_file = just_file.split(".txt")[0] + ".bw"
    bg_file = just_file.split(".txt")[0] + ".bedgraph"
    print(bg_file)
    with open(file_path) as file:
        for line in file:
            cell = line.split("\t")

            score = float(cell[3])
            # print(cell, score)
            ### MAX Score
            if score != 0:
                if cell[0] not in ce_obj:
                    ce_obj[cell[0]] = score
                if cell[1] not in ce_obj:
                    ce_obj[cell[1]] = score
                if score > ce_obj[cell[0]]:
                    ce_obj[cell[0]] = score
                if score > ce_obj[cell[1]]:
                    ce_obj[cell[1]] = score

    out_lines = []
    for ce in ce_obj:
        chrom, c_range = ce.split(":")
        start, end = c_range.split("-")
        stat = ce_obj[ce]
        out_line = chrom + "\t" + start + "\t" + end + "\t" + str(stat) + "\n"
        out_lines.append(out_line)
    print(out_lines)
    with open("ce_score_messy.bedgraph", "w") as file:
        file.writelines(out_lines)
    os.system(f"LC_COLLATE=C sort -k1,1 -k2,2n ce_score_messy.bedgraph > {just_path}/{bg_file}")
    with open(f"{just_path}/{bg_file}") as file:
        for line in file:
            print(line)
    os.system(f"./bedGraphToBigWig {just_path}/{bg_file} mm10Sort.genome {just_path}/{bw_file}")
    # clean-up
    # os.remove("ce_score.bedgraph")
    os.remove("ce_score_messy.bedgraph")


# /Users/blake/Documents/gargLab/for_github/ranking_fisher_1/ranked_pipeline/hg38_sort.chrom
main("/Users/blake/Documents/gargLab/for_github/ranking_fisher_1/ranked_pipeline/data_GM12878/ce_rank.txt")
# main("/Users/blake/Documents/gargLab/fire_story_v3/data/cluster_3_ce_score.txt")
# main("/Users/blake/Documents/gargLab/fire_story_v3/data/cluster_2_ce_score.txt")
# main("/Users/blake/Documents/gargLab/fire_story_v3/data/cluster_1_ce_score.txt")
