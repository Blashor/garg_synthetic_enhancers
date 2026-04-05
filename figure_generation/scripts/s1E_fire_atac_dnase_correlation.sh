#!/bin/bash

#SBATCH --job-name=fire_story_v3
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=scripts/sbatchOutput/firev3_%j.out
#SBATCH --error=scripts/sbatchOutput/firev3_%j.err
#SBATCH --mail-user=blake.torres@yale.edu
#SBATCH --mail-type=ALL

# example script for processing ATACseq
### Environment built with the following:
# mamba create --name gargLab -c conda-forge -c bioconda trim-galore fastx_toolkit fastp bowtie2 samtools deeptools bedtools genrich subread hisat2 r-base 


cd /home/bmt26/scripts/fire_story_v3/



module load miniconda
conda activate gargLab
#python bin128_cluster.py
#python ce_rank_to_metagene_cluster.py 
cd /home/bmt26/garg/FIRE_snakemake/FIRE/results/yaleFiberJul29_2024/trackHub/bw
multiBigwigSummary bins -b all.fire.coverage.bw SRX10040677_dnase.bw SRX23682289_atac.bw -bs 10000 --smartLabels -p max -o dnase_atac_fire_10000.npz

plotCorrelation -in dnase_atac_fire_10000.npz -c spearman -p scatterplot -o plot_10000_s.svg
plotCorrelation -in dnase_atac_fire_10000.npz -c pearson -p scatterplot -o plot_10000_p.svg