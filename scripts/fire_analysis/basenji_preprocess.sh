#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --job-name=basenji_cluster
#SBATCH --gpus=rtx3090:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --output=scripts/sbatchOutput/basenji_%j.out
#SBATCH --error=scripts/sbatchOutput/basenji_%j.err
#SBATCH --mail-user=blake.torres@yale.edu
#SBATCH --mail-type=ALL

module load miniconda

conda activate /gpfs/gibbs/project/garg/bmt26/conda_pkgs/basenji_gpu


cd /home/bmt26/garg/basenji/basenji-master/bin

fasta="../data/GRCm38.p6.genome.fa"
wig_table="/home/bmt26/garg/basenji/basenji-master/data/cluster_bin.txt"
limit_bed="ce.bed"
blacklist="../data/gff_genes_blacklist.bed"
#green
preprocess_data="../data/cluster_aug13"

model_params="/home/bmt26/garg/basenji/basenji-master/models/params_cluster.json"
trained_model="../models/cluster_aug13"

test_output="../output/cluster_aug13"
#green
sat_bed="../sat/sox2_super.bed"
sat_length="600"
#yellow
#sat_bed="../sat/sox_span_super.bed"
sat_length="1000"
sat_bed="/home/bmt26/garg/basenji/basenji-master/sat/sox2_131072.bed"
meme_file="../data/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"

#16384, 8192, 4096
 
### Install Info
#orange for basenji_data there is an issue with chrY that needs to be removed from fasta files 
#orange also for this, it calls two other python files basenji_data (_read and _write), need to put python in front of this
#red for basenji_train need to do this: mamba install -c conda-forge cudatoolkit cudnn
#yellow changed limit behavior to bedtools intersect -u for unique hits that don't fully overlap

### DATA PREP cyan
#-b ${blacklist}
#--crop 16384
cat ${model_params}
python basenji_data.py -s 1 -l 131072 --local -o ${preprocess_data} --limit ${limit_bed} -p 128 -t .1 -v .1 -w 128 ${fasta} ${wig_table}

### TRAIN blue
python basenji_train.py -o ${trained_model} ${model_params} ${preprocess_data}

### TEST purple
python basenji_test.py --ai 0,1,2,3 -o ${test_output} --rc --shifts "1,0,-1" ${model_params} ${trained_model}/model_best.h5 ${preprocess_data}

### SAT ANALYSIS violet

#python basenji_sat_bed.py -f ${fasta} -l ${sat_length} -o ${test_output}/sox2_sat_full --rc -t ${preprocess_data}/targets.txt ${model_params} ${trained_model}/model_best.h5 ${sat_bed}
#python basenji_sat_plot.py --png -l ${sat_length} -o ${test_output}/sox2_sat_full/plots -t ${preprocess_data}/targets.txt ${test_output}/sox2_sat_full/scores.h5

#python /home/bmt26/scripts/basenji/FIRE_bed_to_sat.py

#python /home/bmt26/scripts/basenji/FIRE_bed_predictions.py





#python basenji_motifs.py -m ${meme_file} -p 8 -o ${test_output}/motifs ${model_params} ${trained_model}/model_best.h5 ${preprocess_data}
#python basenji_sat_bed.py -f ${fasta} -l ${sat_length} -o ${test_output}/sox2_sat --rc -t ${preprocess_data}/targets.txt ${model_params} ${trained_model}/model_best.h5 ${sat_bed}
#python basenji_predict_bed.py -f ${fasta} --rc -o ${test_output}/sox2_pred2 -t ${preprocess_data}/targets.txt ${trained_model}/params.json ${trained_model}/model_best.h5 ${sat_bed}

