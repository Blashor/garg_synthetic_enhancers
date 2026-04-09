#!/bin/bash
#This script was run with gastruloid samples derived from Shh PEP-edited lines.


#Trim adapters and filter reads
module load Trim_Galore #(v0.6.7)
module load FastQC #(v0.12.1)

trim_galore --fastqc --illumina --paired read1.fastq read2.fastq --output_dir output_directory

#Generate genome index with STAR
module load STAR #(v2.7.11a)

set -e

GENOME=mm10_no_alt_analysis_set_ENCODE.fasta
GTF=gencode.vM4.annotation.gtf
OUTDIR=star_index

THREADS=8
SJDB_OVERHANG=149

mkdir -p $OUTDIR

STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir $OUTDIR \
     --genomeFastaFiles $GENOME \
     --sjdbGTFfile $GTF \
     --sjdbOverhang $SJDB_OVERHANG
     
#Prepare transcript reference for RSEM
module load RSEM #(v.1.3.3)

GENOME=mm10_no_alt_analysis_set_ENCODE.fasta
GTF=gencode.vM4.annotation.gtf
OUTPREFIX=rsem_index

THREADS=8

mkdir -p $(dirname $OUTPREFIX)

rsem-prepare-reference \
    --gtf $GTF \
    --num-threads $THREADS \
    --star \
    $GENOME \
    $OUTPREFIX


#Run STAR and RSEM with pipeline obtained from ENCODE-DCC/long-rna-seq-pipeline/DAC/STAR_RSEM.sh on GitHub

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory

read1=/directory/read1.fq.gz #gzipped fastq file for read1
read2="" #gzipped fastq file for read1
STARgenomeDir=/directory/star_index
RSEMrefDir=/directory/rsem_index/rsem_index
dataType=str_SE # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=8 # number of threads for STAR
nThreadsRSEM=8 # number of threads for RSEM

# executables

STAR=STAR                             
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig  

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"

case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
$STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=60G

mv Aligned.toTranscriptome.out.bam Tr.bam 

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | sort -S $trBAMsortRAM -T ./ ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S $trBAMsortRAM -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
esac

#'rm' Tr.bam

# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd --calc-ci --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac
  
###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem