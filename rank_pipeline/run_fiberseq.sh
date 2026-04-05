#!/bin/bash

#SBATCH --job-name=fiberSeq_msp
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=128G
#SBATCH --output=scripts/sbatchOutput/fiberSeq_%j.out
#SBATCH --error=scripts/sbatchOutput/fiberSeq_%j.err
#SBATCH --mail-user=blake.torres@yale.edu
#SBATCH --mail-type=ALL
cd /home/bmt26/garg/blake/fiberseq/michele_sep1/hifi_reads

module load miniconda
conda activate gargLab
#cd /home/bmt26/garg/FIRE_snakemake/FIRE/results/GM12878/
#ft fire -vv -- —-extract GM12878.fire.bam > acc.model.results.bed
#python approximate_allAccess.py


mmiRef="/gpfs/gibbs/project/garg/shared/catalogs/minimap2_indexes/GRCm38.primary_assembly.genome.mmi"

List="
m84189_250814_210557_s1.hifi_reads
"
#XSGKI_20210730_S64049_PL100189990-1_D01
#XSGKI_20210222_S64018_PL100164255A-1_D01
#XSGKI_20211122_S64018_PL100189990-1_A01
mkdir -p align_r6 align_fiber_r6 #align_fiber_filter

pbmm2 index /home/bmt26/garg/blake/fiberseq/fasta/r6_insert.fa /home/bmt26/garg/blake/fiberseq/fasta/r6_insert.mmi
mmiRef="/home/bmt26/garg/blake/fiberseq/fasta/r6_insert.mmi"
for f in ${List}; do
	(
   ### Subreads -> CCS w/ kinetics 
	#mkdir -p ${f}
	echo $f
   #### Align CCS

	pbmm2 align ${mmiRef} ${f}.bam align_r6/${f}.bam --sort
    #### Fibertools: Predict m6a + nucleosomes
	ft predict-m6a -t 36 -v -b 1 align/${f}.bam align_fiber_r6/${f}.bam
	bamCoverage -b align_r6/${f}.ccs.bam -o align_r6/${f}.bw
	#ft fire -e -t 36 -v align_fiber_r6/${f}.bam align_fiber_r6/${f}_fire.bed 
	#Extract to bed
	#"""-allow1bpOverlap"""
	ft extract -v -r --msp align_fiber_r6/${f}_msp.bed align_fiber_r6/${f}.bam
	samtools view -H align_fiber_r6/${f}.bam > chroms_samtools.txt
	awk '$1 ~ /^chr[1-9][0-9]*$/ { print }' align_fiber_r6/${f}_msp.bed > align_fiber_r6/${f}_chr.bed
	rm align_fiber_r6/${f}_msp.bed
	#Sort bed
   sortBed -i align_fiber_r6/${f}_chr.bed > align_fiber_r6/${f}_sort.bed
   rm align_fiber_r6/${f}_chr.bed
   #Create bigBed for viewing
   /home/bmt26/scripts/FIBERseq/bedToBigBed/bedToBigBed -allow1bpOverlap align_fiber_r6/${f}_sort.bed /home/bmt26/scripts/FIBERseq/bedToBigBed/mm10.chrom.sizes align_fiber_r6/${f}_msp.bb

	) &

	#
done

wait


