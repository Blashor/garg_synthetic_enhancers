import sys
import os
import re
import json


#
# Generates intergenic whitelist regions
#


def get_whitelist_peaks(gff_file, peaks, whitelist_peaks, chrom_file, promoter_size=500):
    lines = []
    genetypes = []
    out_lines = []
    with open(gff_file) as file:
        for line in file:
            cell = line.split("\t")
            if len(cell) == 9:
                if cell[2] == "gene":
                    geneType = cell[8].split("gene_type=")[1].split(";")[0]
                    geneName = cell[8].split("gene_name=")[1].split(";")[0]
                    if geneType != "lncRNA" and "pseudogene" not in geneType:
                        if cell[6] == "+":
                            cell[3] = int(cell[3]) - promoter_size
                        elif cell[6] == "-":
                            cell[4] = int(cell[4]) + promoter_size
                        if re.search(r"^chr(1\d|\d|X)$", cell[0]):
                            out_lines.append(
                                cell[0] + "\t" + str(cell[3]) + "\t" + str(cell[4]) + "\t" + geneName + "\n"
                            )

    with open("temp/gff_genes.bed", "w") as file:
        file.writelines(out_lines)
    os.system("LC_COLLATE=C sort -k1,1 -k2,2n temp/gff_genes.bed > temp/gff_genes_sort.bed")

    os.system(f"bedtools complement -i temp/gff_genes_sort.bed -g {chrom_file} -L > temp/gff_whitelist.bed")

    os.system(f"bedtools intersect -sorted -F 1 -a temp/gff_whitelist.bed -b {peaks} > {whitelist_peaks}")
    #
    os.remove("temp/gff_genes.bed")
    os.remove("temp/gff_genes_sort.bed")
    os.remove("temp/gff_whitelist.bed")

    return whitelist_peaks


def load_peaks(gff_filtered_file):
    lines = []
    peaks_chr = {}
    with open(gff_filtered_file, "r") as file:
        for line in file:
            cell = line.split("\t")
            if len(cell) == 3:
                if not cell[0] in peaks_chr:
                    peaks_chr[cell[0]] = []
                peaks_chr[cell[0]].append((int(cell[1]), int(cell[2]), "peak"))
    return peaks_chr


def generate_gff_roadblocks(gff3):
    lines = []
    peaks_chr = {}
    with open(gff3, "r") as file:
        lines = file.read().split("\n")
    for line in lines:
        cell = line.split("\t")
        if len(cell) == 9:
            if cell[2] == "gene":
                geneType = cell[8].split("gene_type=")[1].split(";")[0]
                if geneType != "lncRNA" and "pseudogene" not in geneType:
                    if not cell[0] in peaks_chr:
                        peaks_chr[cell[0]] = []
                    peaks_chr[cell[0]].append((int(cell[3]), int(cell[3]), "gff"))
                    peaks_chr[cell[0]].append((int(cell[3]), int(cell[4]), "gff"))
                    peaks_chr[cell[0]].append((int(cell[4]), int(cell[4]), "gff"))
    return peaks_chr


def combine_chr(peak_chr, gff_chr):
    peaks_chr = {}
    for chrom in peak_chr:
        oneList = peak_chr[chrom] + gff_chr[chrom]
        oneList.sort(key=lambda tup: tup[0])
        peaks_chr[chrom] = oneList
    return peaks_chr


def toBedFile(file, stitch_obj):
    bedLines = []
    for stitched in stitch_obj:
        chrom = stitched[0][0]
        start = stitched[0][1]
        end = stitched[-1][2]
        bedLines.append(chrom + "\t" + start + "\t" + end + "\n")
    f = open(file, "w")
    f.writelines(bedLines)
    f.close()


def stitch_peaks(whitelist_peaks, stitch_filename, gff_file, stitch_window=12500):
    # this function generates the integenic filtered bed file used in load_peaks
    peaks_chr = combine_chr(load_peaks(whitelist_peaks), generate_gff_roadblocks(gff_file))
    allStitched = []

    for chrom in peaks_chr:
        latestStitched = []
        for itr in range(len(peaks_chr[chrom]) - 1):
            peak = peaks_chr[chrom][itr]
            if peak[2] == "peak":
                latestStitched.append((chrom, str(peak[0]), str(peak[1])))
                nextPeak = peaks_chr[chrom][itr + 1]
                distance = nextPeak[0] - peak[1]
                # start a newEnhancer
                if distance > stitch_window or nextPeak[2] == "gff":
                    allStitched.append(latestStitched)
                    latestStitched = []
    toBedFile(stitch_filename, allStitched)


def main(gff_file, chrom_file, peaks_file, output):
    # Inputs

    # Generate Intergenic FIRE peaks
    os.system(f"mkdir -p {output}")
    os.system(f"mkdir -p temp")
    whitelist_peaks = get_whitelist_peaks(gff_file, peaks_file, f"{output}/FIRE_peaks_intergenic.bed", chrom_file)
    # Stitch Intergenic FIRE peaks
    stitch_peaks(whitelist_peaks, f"{output}/FIRE_stitched.bed", gff_file)

    return f"{output}/FIRE_stitched", f"{output}/FIRE_peaks_intergenic"


# main()
