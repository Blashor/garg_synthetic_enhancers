import re
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

fire_peaks = "../data/3runs-FDR-FIRE-peaks.bed"
fire_intergenic = "../data/3runV2_intergenic.bed"
fire_genic = "../data/3runV2_genic.bed"
gff = "../data/gencode.vM23.annotation.gff3"
os.system(f"bedtools merge -i {fire_peaks} > {fire_peaks}.merge.bed")


intergenic_peaks = set()

with open(fire_intergenic) as file:
    for line in file:
        cell = line.strip().split("\t")
        chrom = f"{cell[0]}:{cell[1]}-{cell[2]}"
        intergenic_peaks.add(chrom)
print(len(intergenic_peaks))

lines = []
with open(f"{fire_peaks}.merge.bed") as file:
    for line in file:
        cell = line.strip().split("\t")
        chrom = f"{cell[0]}:{cell[1]}-{cell[2]}"
        chrom_line = f"{cell[0]}\t{cell[1]}\t{cell[2]}"
        if chrom not in intergenic_peaks:
            lines.append(chrom_line + "\n")

with open(fire_genic, "w") as file:
    file.writelines(lines)
itr = 0


def gff_to_bed(gff="../data/gencode.vM23.annotation.gff3"):
    gff_bed_lines = []
    with open(gff) as file:
        for line in file:
            cell = line.split("\t")
            if "lncRNA;" in line:
                gff_id = cell[2] + "_lncRNA"
                gff_bed_lines.append(f"{cell[0]}\t{cell[3]}\t{cell[4]}\t{gff_id}\n")
            elif "pseudogene;" not in line and len(cell) > 2:
                gff_id = cell[2] + "_protein_coding"
                gff_bed_lines.append(f"{cell[0]}\t{cell[3]}\t{cell[4]}\t{gff_id}\n")
                if "gene" in cell[2]:
                    if cell[6] == "+":
                        gff_bed_lines.append(
                            f"{cell[0]}\t{max(0,int(cell[3])-500)}\t{cell[3]}\tpromoter_protein_coding\n"
                        )
                        gff_bed_lines.append(
                            f"{cell[0]}\t{max(0,int(cell[3])-5000)}\t{cell[3]}\t5kb_proximal_protein_coding\n"
                        )
                    elif cell[6] == "-":
                        gff_bed_lines.append(f"{cell[0]}\t{cell[4]}\t{int(cell[4])+500}\tpromoter_protein_coding\n")
                        gff_bed_lines.append(
                            f"{cell[0]}\t{cell[4]}\t{int(cell[4])+5000}\t5kb_proximal_protein_coding\n"
                        )

            else:
                pass
    with open("../data/vm23_gff.bed", "w") as file:
        file.writelines(gff_bed_lines)
    lines = []
    linesToFile = []
    with open("../../../CHIPseq/H3K27ac.narrowpeak", "r") as file:
        lines = file.read().split("\n")
    for line in lines:
        cells = line.split("\t")[:3]
        if cells[0][:3] == "chr":
            cells.append("h3k27ac_peak")
            linesToFile.append("\t".join(cells) + "\n")
    gffBed = open("../data/vm23_gff.bed", "a")
    gffBed.writelines(linesToFile)
    gffBed.close()
    os.system("LC_COLLATE=C sort -k1,1 -k2,2n ../data/vm23_gff.bed > ../data/vm23_gff.sort.bed")
    os.remove("../data/vm23_gff.bed")

    os.system(
        f"bedtools intersect -sorted -a {fire_intergenic} -b ../data/vm23_gff.sort.bed  -wao > ../data/fire_interGFF.bed"
    )
    os.system(
        f"bedtools intersect -sorted -a {fire_genic} -b ../data/vm23_gff.sort.bed  -wao > ../data/fire_geneGFF.bed"
    )


allowed = set(
    [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chr23",
        "chrX",
    ]
)


def load_FIREgff():
    lines = []
    fireOverlap = {}
    for file_name in ["../data/fire_interGFF.bed", "../data/fire_geneGFF.bed"]:
        with open(file_name, "r") as file:
            lines = file.read().split("\n")
        for line in lines:
            cells = line.split("\t")
            if cells[0] in allowed:
                if len(cells) > 6:
                    fireId = cells[0] + ":" + cells[1] + "-" + cells[2]
                    if fireId in fireOverlap:
                        pass
                    else:
                        fireOverlap[fireId] = []
                    # print(cells[6])
                    if file_name == "../data/fire_geneGFF.bed":
                        if cells[6] in hierarchy_genic_allowed:
                            fireOverlap[fireId].append(cells[6])
                    else:
                        fireOverlap[fireId].append(cells[6])
                else:
                    pass
                    print(cells)
    return fireOverlap
    # fireOverlap[]


hierarchy = {
    "gene_lncRNA": 0,
    "five_prime_UTR_protein_coding": 1,
    "three_prime_UTR_protein_coding": 2,
    "exon_protein_coding": 3,
    "gene_protein_coding": 4,
    "promoter_protein_coding": 5,
    "5kb_proximal_protein_coding": 6,
    "h3k27ac_peak": 7,
    ".": 8,
}
hierarchy_genic_allowed = {
    "five_prime_UTR_protein_coding": 1,
    "three_prime_UTR_protein_coding": 2,
    "exon_protein_coding": 3,
    "gene_protein_coding": 4,
    "promoter_protein_coding": 5,
}
label_groups = {
    1: [],
    2: [],
    3: [],
    4: [],
    5: [],
    6: [],
    8: [],
    7: [],
    0: [],
}
labels = [
    "5' UTR",
    "3' UTR",
    "Exons",
    "Introns",
    "Promoters",
    "5kb Proximal",
    "Intergenic",
    "H3K27ac",
    "lncRNA",
]


def load_FIRE_peaks():
    lines = []
    with open("../data/fire_interGFF.bed", "r") as file:
        lines = file.read().split("\n")
    return len(lines) - 1


def main():
    plt.figure()
    # groupData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    groupData = []
    sns.set()
    sns.set_palette("gray")  # Set grayscale color palette
    fireOverlap = load_FIREgff()
    for p_i, peak in enumerate(fireOverlap):
        priority = 8
        for peakInfo in fireOverlap[peak]:
            if peakInfo in hierarchy:
                hierarchyNum = hierarchy[peakInfo]
                if hierarchyNum < priority:
                    priority = hierarchyNum
        label_groups[priority].append(peak)
        # print(p_i, priority, peak)
    totalFire = load_FIRE_peaks()

    for arr in label_groups:
        groupData.append(len(label_groups[arr]))
    print(groupData)
    itr = 0
    outlines = []
    for arr in label_groups:
        for p_i, p in enumerate(label_groups[arr]):
            chrom, coords = p.split(":")
            start, end = coords.split("-")
            outlines.append(f"{chrom}\t{start}\t{end}\t{labels[arr]}\n")
            # print(itr, p, labels[arr])
            # itr += 1
    with open("gff_distribution.txt", "w") as file:
        file.writelines(outlines)
    total_sum = sum(groupData)
    percentages = [value / total_sum * 100 for value in groupData]
    vals = [value for value in groupData]
    combined_labels = [f"{label}:\n{val:,} ({percent:.1f}%)" for label, percent, val in zip(labels, percentages, vals)]

    wedges, texts = plt.pie(
        groupData,
        labels=combined_labels,  # Use the modified labels with percentages
        autopct=None,  # Disable autotext since percentages are in labels
        pctdistance=1.15,
        labeldistance=1.36,
        startangle=0,
        textprops=dict(ha="center"),
    )

    groups = [[0, 1, 2, 3, 4], [5, 6, 7, 8]]

    radfraction = 0.1
    for group in groups:
        ang = np.deg2rad((wedges[group[-1]].theta2 + wedges[group[0]].theta1) / 2)
        for j in group:
            center = radfraction * wedges[j].r * np.array([np.cos(ang), np.sin(ang)])
            wedges[j].set_center(center)
            texts[j].set_position(np.array(texts[j].get_position()) + center)
            # percs[j].set_position(np.array(percs[j].get_position()) + center)
    # ax.autoscale(True)
    plt.show()
    plt.savefig("../figures/s2b_pie.svg")
    plt.show()


gff_to_bed()
main()
