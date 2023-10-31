#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -N 1
#SBATCH -t 24:00:00
#SBATCH -J CLUMPP
#SBATCH -M rackham

ml R/3.6.1
ml R_packages/3.6.1

Rscript ../scripts/Summarise_AdmixureRuns.R
