#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --job-name=INSIDE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=rtx5000:1
#SBATCH --time=2-00:00:00
#SBATCH --mem=40G
#SBATCH --output=scripts/sbatchOutput/INSIDE_%j.out
#SBATCH --error=scripts/sbatchOutput/INSIDE_%j.err
#SBATCH --mail-type=ALL



module load miniconda

conda activate /gpfs/gibbs/project/garg/bmt26/conda_pkgs/basenji_gpu

cd /home/bmt26/garg/basenji/basenji-master/bin


 
### Install Info

#orange for basenji_data there is an issue with chrY that needs to be removed from fasta files 
#orange also for this, it calls two other python files basenji_data (_read and _write), need to put python in front of this
#red for basenji_train need to do this: mamba install -c conda-forge cudatoolkit cudnn
#yellow changed limit behavior to bedtools intersect -u for unique hits that don't fully overlap

#json setup files should be placed in /home/bmt26/garg/basenji/basenji-master/json

###

python /home/bmt26/scripts/INSIDE/INSIDE.py cluster140
#python /home/bmt26/scripts/INSIDE/INSIDE_score_parents_no_selection.py cluster140
#python /home/bmt26/scripts/INSIDE/INSIDE_within_genome.py 1000o14 
