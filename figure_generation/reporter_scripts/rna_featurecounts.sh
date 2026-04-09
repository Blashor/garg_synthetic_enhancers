#!/bin/bash

#Trim adapters and filter reads
module load Trim_Galore #(v0.6.7)
module load FastQC #(v0.12.1)

trim_galore --fastqc --illumina --paired read1.fastq read2.fastq --output_dir output_directory

#Generate genome index with STAR
module load STAR #(v2.7.11a)

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /directory/star_index \
  --genomeFastaFiles GRCm38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.vm23.grcm38.p6.primary.anno.gtf \
  --sjdbOverhang 149

#Align processed paired-end reads to STAR genome index
STAR \
  --runThreadN 4 \
  --genomeDir /directory/star_index \
  --readFilesIn /directory/read1_processed.fq /directory/read2_processed.fq \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix /directory/sample_align
  
#Count mapped reads with featureCounts
module load Subread #(v2.0.3)

featureCounts -T 8 \
    -a gencode.vm23.grcm38.p6.primary.anno.gtf \
    -p -B -C \
    -o counts.txt \
    sample_alignAligned.out.bam
    
cat counts.txt | tr '\t' ',' > counts.csv
    
#Convert from Mouse (Mus musculus) ENSEMBL Gene IDs to gene symbols in counts.txt
#1. Generate a list with ENSEMBL Gene IDs and corresponding gene symbols
awk '$3 == "gene" { 
    match($0, /gene_id "[^"]+"/, gid); 
    match($0, /gene_name "[^"]+"/, gname); 
    if (gid[0] && gname[0]) { 
        gsub(/gene_id "|"/, "", gid[0]); 
        gsub(/gene_name "|"/, "", gname[0]); 
        print gid[0] "\t" gname[0] 
    } 
}' gencode.vm23.grcm38.p6.primary.anno.gtf > gene.tsv

#2. Switch ENSEMBL Gene IDs to gene names
awk -F',' -v OFS=',' '
BEGIN {
    while ((getline < "gene.tsv") > 0) {
        split($0, a, "\t")
        map[a[1]] = a[2]
    }
}
NR == 1 {
    print
    next
}
{
    if ($1 in map) {
        $1 = map[$1]
    }
    print
}
' counts.csv > counts_genesymbol.csv