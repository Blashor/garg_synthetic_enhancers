#!/bin/bash

#SBATCH --partition=ycga
#SBATCH --job-name=run_FIRE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=scripts/sbatchOutput/run_FIRE_%j.out
#SBATCH --error=scripts/sbatchOutput/run_FIRE_%j.err
#SBATCH --mail-user=blake.torres@yale.edu
#SBATCH --mail-type=ALL


module load miniconda
conda activate snakemake
#export PATH=$PATH:/home/bmt26/garg/FIRE_snakemake/ucsc
#export SNAKEMAKE_CONDA_PREFIX=/gpfs/gibbs/project/garg/bmt26/conda_envs/snakemake
#samtools index /home/bmt26/garg/fiberSeq/yaleFiberJun17_2024/align_fiber/m84189_240614_003708_s4.bam


#MERGING
#cd /home/bmt26/garg/fiberSeq/Jul29/align_fiber

#one=XSGKI_20210222_S64018_PL100164255A-1_D01.bam
#two=XSGKI_20210730_S64049_PL100189990-1_D01.bam
#three=XSGKI_20211122_S64018_PL100189990-1_A01.bam

#samtools index ${one}
#samtools index ${two}
#samtools index ${three}
#samtools merge merged_unsort.bam ${one} ${two} ${three}
#samtools sort -o merge3.bam merged_unsort.bam
#samtools index merge3.bam

#rm merged_unsort.bam

cd /home/bmt26/garg/FIRE_snakemake/FIRE
./fire --configfile /home/bmt26/garg/FIRE_snakemake/FIRE/config/config.yaml --unlock
./fire --configfile config/config.yaml