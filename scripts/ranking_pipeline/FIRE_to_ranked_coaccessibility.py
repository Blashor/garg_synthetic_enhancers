import os
import sys
import argparse
import FIRE_stitcher
import objPrep


#
# FIRE Ranking Pipeline
# Outputs FIRE peak pairs ranked by co-accessibility
#
def build_pre_matrix(stitch_file, peak_file, acc_fiber_bed, output_name):
    # stitch_file = "/home/bmt26/garg/fiberSeq/FIRE/SE/FIRE_peaks_stitch_ROSElike"
    # output_name = "FIRE_peaks_stitch_ROSElike"
    # peak_file = "/home/bmt26/garg/fiberSeq/FIRE/SEconstituents/FIRE_peaks"
    # acc_fiber_bed = "/home/bmt26/garg/fiberSeq/FIRE/FIRE_fiber_calls/acc.model.results.bed"

    # Sort
    # print(f"{stitch_file}.bed")
    os.system(f"LC_COLLATE=C sort -k1,1 -k2,2n {stitch_file}.bed > {stitch_file}_sort.bed")
    os.system(f"LC_COLLATE=C sort -k1,1 -k2,2n {peak_file}.bed > {peak_file}_sort.bed")

    # Intersect stitch + enh
    os.system(
        f"bedtools intersect -sorted -F 1 -a {stitch_file}_sort.bed -b {peak_file}_sort.bed -wa -wb | awk -F'\\t' '{{OFS=\"\\t\"; print $4, $5, $6, $1, $2, $3}}' > {stitch_file}_overlap.bed"
    )

    # Intersect stitch + enh + acc.bed
    os.system(
        f"bedtools intersect -sorted -a {stitch_file}_overlap.bed -b {acc_fiber_bed} -wo > {output_name}/preAwk_Cov.bed"
    )
    os.system(
        f"awk -F'\\t' '{{OFS=\"\\t\"; print $1,$2,$3,$4,$5,$6,$10,$16,$NF}}' {output_name}/preAwk_Cov.bed > {output_name}/Cov.bed"
    )

    # Creates object representation of above data
    # os.system(f"python3 objPrep.py ${output_name}")

    os.system(f"rm {stitch_file}_sort.bed")
    os.system(f"rm {peak_file}_sort.bed")
    os.system(f"rm {stitch_file}_overlap.bed")
    os.system(f"rm {output_name}/preAwk_Cov.bed")


def main():
    parser = argparse.ArgumentParser(description="FIRE output to ranked co-accessibility")

    parser.add_argument("-p", "--peak", type=str, help="Path to the peaks file", required=True)
    parser.add_argument("-c", "--chrom", type=str, help="Path to the chromosome sizes file", required=True)
    parser.add_argument(
        "-a", "--acc", type=str, help="Path to SORTED acc.model.results.bed (FIRE output)", required=True
    )
    parser.add_argument("-g", "--gff", type=str, help="Path to the gff3 file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Provide an output name (String)", required=True)

    args = parser.parse_args()

    # FIRE stitcher
    stitched_name, intergenic_peaks = FIRE_stitcher.main(args.gff, args.chrom, args.peak, args.output)
    print(stitched_name, intergenic_peaks)

    # Build object
    build_pre_matrix(stitched_name, intergenic_peaks, args.acc, args.output)
    objPrep.main(args.output)
    # Rank object
    os.system("mkdir -p data")
    os.system(f"python3 fireAB_ranking.py {args.output}")


# acc_file="/home/bmt26/garg/Michele/fiberSeq/Nov7/fire_snakemake/results/yaleFiberNov7_2024/fire/acc.model.results.sort.bed"
"""
acc_dir="/home/bmt26/garg/blake/FIRE_ranking/jun27/"
k562="/home/bmt26/garg/blake/FIRE_ranking/3runs_v2/"
acc_file="/home/bmt26/garg/FIRE_snakemake/FIRE/results/sep1/fiber-calls/approxAccess.sort.bed"
peaks=${k562}3runs_2_merge.bed
gff=${k562}gencode.vM23.annotation.gff3
chrom_sizes=${k562}GRCm38.primary_assembly.sort.chrom
output_name="sep1_run"
python FIRE_to_ranked_coaccessibility.py -p ${peaks} -g ${gff} -a ${acc_file} -o ${output_name} -c ${chrom_sizes}
"""
main()
