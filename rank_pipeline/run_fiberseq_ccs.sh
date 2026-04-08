#!/bin/bash

#SBATCH --job-name=fiberSeq_msp
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=scripts/sbatchOutput/fiberSeq_%j.out
#SBATCH --error=scripts/sbatchOutput/fiberSeq_%j.err
#SBATCH --mail-user=blake.torres@yale.edu
#SBATCH --mail-type=ALL
cd /home/bmt26/garg/fiberSeq


sifExec="apptainer exec /home/bmt26/scripts/FIBERseq/fiberseq.sif"
sifExecGPU="apptainer exec --nv /home/bmt26/scripts/FIBERseq/fiberseq.sif"
mmiRef="/gpfs/gibbs/project/garg/shared/catalogs/minimap2_indexes/GRCm38.p6.genome.mmi"
List="
XSGKI_20210730_S64049_PL100189990-1_D01
XSGKI_20211122_S64018_PL100189990-1_A01
XSGKI_20210222_S64018_PL100164255A-1_D01
"

#XSGKI_20210730_S64049_PL100189990-1_D01
#XSGKI_20211122_S64018_PL100189990-1_A01

#Still fiberseqing
#XSGKI_20210222_S64018_PL100164255A-1_D01
echo "this one"
mkdir -p align align_fiber #align_fiber_filter
#
chunkNum=24
for f in ${List}; do
	(
   ### Subreads -> CCS w/ kinetics 
	mkdir -p ${f}
	${sifExec} pbindex ${f}.subreads.bam
##
	for ((i=1; i<=${chunkNum}; i++)); do
	    chunk="${i}/${chunkNum}"
	    echo "${chunk}"
	   ${sifExec} ccs ${f}.subreads.bam ${f}/${f}.ccs.${i}.bam --chunk ${chunk} --hifi-kinetics --log-level INFO --metrics-json ${f}/${f}.${i}_metrics.txt --report-file ${f}/${f}.${i}_report.txt &
   #done
	wait
	${sifExec} pbmerge -o ${f}.ccs.bam ${f}/${f}.ccs.*.bam
	${sifExec} pbindex ${f}.ccs.bam
#
   #### Align CCS
	${sifExec} pbmm2 align ${mmiRef} ${f}.ccs.bam align/${f}.ccs.bam --sort
   #### Fibertools: Predict m6a + nucleosomes
	${sifExec} ft predict-m6a -v -b 1 align/${f}.ccs.bam align_fiber/${f}.bam
	#red
	#Bug workaround for our current version of fiberseq, believed to be fixed in current version vvvvvvv
	${sifExec} ft add-nucleosomes -v align_fiber/${f}_bug.bam align_fiber/${f}.bam
	#Extract to bed
	${sifExec} ft extract -v -r --msp align_fiber/${f}_msp.bed align_fiber/${f}.bam
	awk '$1 ~ /^chr[1-9][0-9]*$/ { print }' align_fiber/${f}_msp.bed > align_fiber/${f}_chr.bed
	rm align_fiber/${f}_msp.bed
	#Sort bed
    LC_COLLATE=C sort -k1,1 -k2,2n align_fiber/${f}_chr.bed > align_fiber/${f}_msp2.bed
    rm align_fiber/${f}_chr.bed
    #Create bigBed for viewing
    /home/bmt26/scripts/FIBERseq/bedToBigBed/bedToBigBed align_fiber/${f}_msp2.bed /home/bmt26/scripts/FIBERseq/bedToBigBed/mm10.chrom.sizes align_fiber/${f}_msp2.bb

    # do other visualizations, less than 5kb cut

    #samtools index align_fiber/${f}.bam
	#samtools view -b align_fiber/${f}.bam -q 30 > align_fiber_filter/${f}_filtered.bam
   	#samtools sort align_fiber_filter/${f}_filtered.bam -o align_fiber_filter/${f}_coord.bam
   	#samtools index align_fiber_filter/${f}_coord.bam
    #${sifExec} ft extract -v --m6a align_fiber_filter/${f}_met.bed align_fiber_filter/${f}_coord.bam
    #bamCoverage -b align_fiber_filter/${f}_coord.bam -o align_fiber_filter/${f}.bw --binSize 5 --numberOfProcessors max/2
	) &

	#
done

wait


