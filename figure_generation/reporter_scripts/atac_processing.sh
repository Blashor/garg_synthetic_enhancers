#!/bin/bash

#Trim adapters and filter reads
module load Trim_Galore #(v0.6.7)
module load FastQC #(v0.12.1)

trim_galore --cores 4 --fastqc --nextera --paired read1.fastq read2.fastq --output_dir outputdirectory

#Generation of custom genomes for reporter-derived data
#Append fasta file of enhancer pair reporter sequence to that of reference genome
cat reporter.fa GRCm38.primary_assembly.genome.fa > custom_genome.fa

#Generation of custom genomes for Meox1 and Shh PEP-derived data
#1. Extract chromosome 11 (Meox1) or 5 (Shh) from reference genome
input="GRCm38.primary_assembly.genome.fa"
output="chr11.fa"
chromosome="chr11 11" #input chromosome of interest

awk -v chr=">$chromosome" '
    $0 ~ /^>/ {printing = ($0 == chr)}
    printing {print}
' "$input" > "$output"

#2. Rewrite chromosome 11 (Meox1) or 5 (Shh) in text editor to replace endogenous sequence with PEP sequences
#3. Swap the edited chromosome with the unedited chromosome in the reference genome
module load miniconda
conda activate conda_environment #(-c conda-forge -c bioconda trim-galore fastx_toolkit fastp bowtie2 samtools deeptools bedtools genrich subread hisat2 r-base)

cd /directory/to/script
python chromosome_switch.py

#Index custom genomes with Bowtie2
module load SAMtools #(v1.21)
module load Bowtie2 #(v2.5.1)

samtools faidx custom_genome.fa
bowtie2-build custom_genome.fa customindex

#Align processed paired-end data to custom genomes with Bowtie2
bowtie2 -p 8 --very-sensitive -k 10 -x customindex -q -1 read1_processed.fq -2 read2_processed.fq -S sample.sam

#Use SAMtools to convert the SAM file to a BAM file and subsequently sort and index it
samtools view -S -b sample.sam > sample.bam
samtools sort sample.bam -o sample.sort.bam
samtools index sample.sort.bam

#Use deepTools to normalize the sorted and indexed BAM file and generate a bigWig file for visualization
module load deepTools #(v3.5.5)

bamCoverage --normalizeUsing CPM --binSize 10 -b sample.sort.bam -o sample.bw 

#Visualize all reporter sequences in one view in IGV
module load miniconda
conda activate conda_environment #(-c conda-forge -c bioconda trim-galore fastx_toolkit fastp bowtie2 samtools deeptools bedtools genrich subread hisat2 r-base)

cd /directory/to/script
python oneview.py