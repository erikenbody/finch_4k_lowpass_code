#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J Admixture_Projection
#SBATCH -M rackham

#################################################################################################################
# PERFORM ADMIXTURE ANALYSIS WITH PROJECTION
# A. Sendell-price, May 2022

# This script can be run as follows:
# sbatch scripts/RunAdmixtureProjection.sh VCF_PATH TRAINING_SAMPLE_INFO PROJECT_SAMPLE_INFO OUT_PREFIX K_MIN K_MAX SCRIPTS_DIR
# VCF_PATH = path to VCF file containing all high coverage and low coverage phased genotypes
# TRAINING_SAMPLE_INFO = list of sample IDs to be used for the initial admixture run (no. per spp. roughly equal)
# PROJECT_SAMPLE_INFO = list of sample IDs to be used for the projection admixture runs (all samples of interest)
# OUT_PREFIX = prefix to be used for outputs
# K_VAL = K value to test
# SCRIPTS_DIR = path to directory containing "Admixture_Launcher.sh" script

#################################################################################################################

#Load conda environment and modules
module load bioinfo-tools bcftools vcftools plink ADMIXTURE/1.3.0

#Get variable names from sbatch submission
VCF_PATH=${1}
TRAINING_SAMPLE_INFO=${2}
PROJECT_SAMPLE_INFO=${3}
OUT_PREFIX=${4}
K_MIN=${5}
K_MAX=${6}
SCRIPTS_DIR=${7}

#Create output directory and move into it
mkdir $OUT_PREFIX
cd $OUT_PREFIX

#Save variables to file for record
date > ${OUT_PREFIX}.log
echo " " >> ${OUT_PREFIX}.log
echo "Variables passed were as follows:" >> ${OUT_PREFIX}.log
echo "VCF: " $VCF_PATH >> ${OUT_PREFIX}.log
echo "Training samples: " $TRAINING_SAMPLE_INFO >> ${OUT_PREFIX}.log
echo "Projected samples: " $PROJECT_SAMPLE_INFO >> ${OUT_PREFIX}.log
echo "K_min " $K_MIN >> ${OUT_PREFIX}.log
echo "K_max " $K_MAX >> ${OUT_PREFIX}.log

#Create sample lists from info files
cat $TRAINING_SAMPLE_INFO | cut -f 1 > trainAdmix.list
cat $PROJECT_SAMPLE_INFO | cut -f 1 > projAdmix.list
cp $PROJECT_SAMPLE_INFO sample_info.tmp

#Create directory for input files
mkdir admixture_input

#Subset VCF keeping only samples in lists. Create two versions, one containing training samples
#the other containing all focal samples
vcftools --gzvcf $VCF_PATH --keep trainAdmix.list --stdout --recode \
| bgzip > admixture_input/trainAdmix.vcf.gz
vcftools --gzvcf $VCF_PATH --keep projAdmix.list --stdout --recode \
| bgzip > admixture_input/projAdmix.vcf.gz

#Extract sample order from VCF file
bcftools query -l admixture_input/projAdmix.recode.vcf.gz > sample_order.tmp

#Convert VCFs to plink bed, bim, fam required by admixture
plink --vcf admixture_input/trainAdmix.vcf.gz --const-fid --autosome-num 30 --make-bed --out admixture_input/trainAdmix --allow-extra-chr
awk '{$1=0;print $0}' admixture_input/trainAdmix.bim >  admixture_input/trainAdmix.bim.tmp
mv admixture_input/trainAdmix.bim.tmp admixture_input/trainAdmix.bim
plink --vcf admixture_input/projAdmix.vcf.gz --const-fid --autosome-num 30 --make-bed --out admixture_input/projAdmix --allow-extra-chr
awk '{$1=0;print $0}'  admixture_input/projAdmix.bim >  admixture_input/projAdmix.bim.tmp
mv admixture_input/projAdmix.bim.tmp admixture_input/projAdmix.bim

#Remove VCF files (these aren't needed anymore)
rm admixture_input/trainAdmix.vcf.gz
rm admixture_input/projAdmix.recode.vcf.gz

#For each K value in range K_MIN:K_MAX do the following
for K_VAL in $(seq $K_MIN $K_MAX)
do
    #Make directory for K value and move into it
    mkdir K${K_VAL}
    cd K${K_VAL}
  
    #Conduct initial admixture run using training sample dataset 
    admixture --cv ../admixture_input/trainAdmix.bed $K_VAL > trainAdmix_K$K_VAL.log

    #Using allele frequencies learnt in the initial run launch 100 independent admixture runs
    #this time using the full sample set. This is done using the script "Admixture_100Run_Launcher.sh"
    sbatch ${SCRIPTS_DIR}Admixture_100Run_Launcher.sh $K_VAL

    cd ../

done



