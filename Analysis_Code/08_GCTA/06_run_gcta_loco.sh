#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 9
#SBATCH -t 3-00:00:00
#SBATCH -J gcta
#SBATCH -e gcta_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gcta_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#NOTE THIS CAN BE MEMORY INTENSE
#even tho multithreading doesnt seem to work still need lots of cores

ml bioinfo-tools qctool/2.0.6 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/37_GCTA
cd $WORK_D

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#############################################################

SUBSET=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_samples.txt
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_gcta.txt
COV=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_weight_covariate_gcta.txt

mkdir -p results_${OUTNAME}

##run with all chromosome GRM
gcta64 --mlma-loco --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --mbfile ${OUTNAME}/path_to_all_bfiles.txt --pheno $PHENO --qcovar $COV --out results_${OUTNAME}/geno_assoc_loco --thread-num 5 --autosome-num 49

#other phenos

#run for weight and PC2
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_PC2_gcta.txt
gcta64 --mlma-loco --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --mbfile ${OUTNAME}/path_to_all_bfiles.txt --pheno $PHENO --qcovar $COV --out results_${OUTNAME}/billPC2_geno_assoc --thread-num 9 --autosome-num 49
#run with weight as covariate and INTs
COV=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_weightINT_covariate_gcta.txt
gcta64 --mlma-loco --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --mbfile ${OUTNAME}/path_to_all_bfiles.txt --pheno $PHENO --qcovar $COV --out results_${OUTNAME}/billPC2_geno_assoc_weightINTcov --thread-num 5 --autosome-num 49

##for weight need to make a new covariates file
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_weightINT_gcta.txt
COV=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_covariate_gcta.txt
gcta64 --mlma-loco --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --mbfile ${OUTNAME}/path_to_all_bfiles.txt --pheno $PHENO --qcovar $COV --out results_${OUTNAME}/weightINT_geno_assoc --thread-num 9 --autosome-num 49
