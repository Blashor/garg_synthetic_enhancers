#!/bin/bash

#The bam files inputted in this script were generated from the script rna_featurecounts.sh.

#Sort and index the bam file with SAMtools and then use deepTools to create a bigWig file for visualization
module load SAMtools #(v1.21)
module load deepTools #(v3.5.5)

samtools sort sample_alignAligned.out.bam -o sample.sort.bam
samtools index sample.sort.bam

bamCoverage --normalizeUsing BPM --binSize 10 -b sample.sort.bam -o sample.bw 