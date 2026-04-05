def main():
    # obj_to_bed()
    # OSN: /Users/blake/Documents/gargLab/FIBERseq/fiberSeqBeds/Analysis/youngSEconstiuentmm10.bed
    # FIRE: FIRE_ce_from_stich_v2.bed
    matName = "3runs_2_merge.gz"
    heatFile = "fire_3runs_BOX.svg"
    beds = "FIRE_stich_v2.bed FIRE_ce_from_stich_v2.bed OSN_q10_young_CHPATLS.bed ../data/Oct_chpAtl.bed ../data/Sox_chpAtl.bed ../data/Nanog_chpAtl.bed OSN_stitch.bed"
    beds = "../data/3runs_2_merge.bed"
    # h3k27ac,ctcf,H3K4me1, med1,p300,
    bws = "../data/SRX027332_H3K27ac.bw ../data/SRX000540_CTCF.bw ../data/SRX000583_H3K4me1.bw ../data/SRX022695_Med1.bw ../data/SRX143842_p300.bw"
    # scale-regions
    deepCompute = f"computeMatrix reference-point --referencePoint center -p 8 -S {bws} -R {beds} -a 1500 -b 1500 --binSize 50 --missingDataAsZero -o {matName}"
    # deepCompute = f"computeMatrix scale-regions -p 8 -S {bws} -R {beds} -a 1500 -b 1500 --binSize 10 --missingDataAsZero -o {matName}"

    # os.system(deepCompute)
    # os.system(f"plotHeatmap -m {matName} -out {heatFile}")
    print(heatFile)
    os.system(f"plotHeatmap --kmeans 3 -m {matName} -out {heatFile}.kmeans3.svg --outFileSortedRegions {heatFile}3.bed")
    # os.system(f"plotHeatmap --kmeans 4 -m {matName} -out {heatFile}.kmeans4.svg --outFileSortedRegions {heatFile}4.bed")


# matrix_clustering()
#
main()
