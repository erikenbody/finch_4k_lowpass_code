#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -N 1
#SBATCH --array=1-100:1
#SBATCH -t 24:00:00
#SBATCH -J Admixture_launcher
#SBATCH -M snowy

#Load conda environment and modules
module load bioinfo-tools ADMIXTURE/1.3.0

#Set K value parameter using variable names from sbatch submission
K_VAL=${1}

#Set run number using slurm array job ID
RUN=$SLURM_ARRAY_TASK_ID

#Make directory for run and move into it
mkdir Run_${RUN}
cd Run_${RUN}

#Create symbolic (soft) links to required files
ln -s ../trainAdmix.${K_VAL}.P projAdmix.run${RUN}.${K_VAL}.P.in
ln -s ../../admixture_input/projAdmix.bed projAdmix.run${RUN}.bed
ln -s ../../admixture_input/projAdmix.bim projAdmix.run${RUN}.bim
ln -s ../../admixture_input/projAdmix.fam projAdmix.run${RUN}.fam

#Run admixture projection analysis with cross-validation mode on
admixture -P projAdmix.run${RUN}.bed $K_VAL --cv > projAdmix.run${RUN}.${K_VAL}.log

