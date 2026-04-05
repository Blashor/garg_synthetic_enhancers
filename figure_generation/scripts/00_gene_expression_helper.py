import os
import json
import sys
import math


def constituent_enhancers_to_bed():
    # does the range from first cluster start to end
    out_lines = []
    with open("../data/3runs_ce_rank.txt") as file:
        for line in file:
            cell = line.split("\t")
            chrom, c_range = cell[0].split(":")
            start, end = c_range.split("-")
            out_lines.append("\t".join((chrom, start, end, cell[4])) + "\n")
            chrom, c_range = cell[1].split(":")
            start, end = c_range.split("-")
            out_lines.append("\t".join((chrom, start, end, cell[4])) + "\n")
    print(out_lines)
    with open("../data/ce_stitch_pre.bed", "w") as file:
        file.writelines(out_lines)
    os.system("LC_COLLATE=C sort -k1,1 -k2,2n ../data/ce_stitch_pre.bed > ../data/ce_stitch.bed")
    os.remove("../data/ce_stitch_pre.bed")


def gtf_gene_to_bed():
    out_lines = []
    with open("../data/gencode.vM23.annotation.gff3") as file:
        for line in file:
            cell = line.split("\t")
            if len(cell) == 9:
                # print(cell[2])
                if cell[2] == "gene":
                    gene_info = {}
                    for pair in cell[8].split(";"):
                        key, value = pair.split("=")
                        gene_info[key] = value
                    if gene_info["gene_type"] == "protein_coding" or "miRNA" in gene_info["gene_type"]:
                        # if gene_info["gene_type"] != "null":
                        if cell[6] == "+":
                            start = cell[3]
                            end = str(int(cell[3]) + 3)
                        else:
                            start = str(int(cell[4]) - 3)
                            end = cell[4]
                        out_lines.append(
                            "\t".join(
                                (
                                    cell[0],
                                    start,
                                    end,
                                    gene_info["gene_id"],
                                    gene_info["gene_name"],
                                )
                            )
                            + "\n"
                        )
    with open("../data/genes_pre.bed", "w") as file:
        file.writelines(out_lines)
    os.system("LC_COLLATE=C sort -k1,1 -k2,2n ../data/genes_pre.bed > ../data/genes_gtf.bed")
    os.remove("../data/genes_pre.bed")


def bedtools_intersect(a, b, c):
    os.system("bedtools intersect -sorted -a " + a + " -b " + b + " -wo > " + c)


def bedtools_closest(a, b, c):
    os.system("bedtools closest -d -a " + a + " -b " + b + " > " + c)
    os.remove(a)
    os.remove(b)


def main():
    gtf_gene_to_bed()
    constituent_enhancers_to_bed()
    bedtools_closest(
        "../data/ce_stitch.bed",
        "../data/genes_gtf.bed",
        "../data/ce_gene_closest.bed",
    )


main()


# clustered_enhancers_to_bed()
