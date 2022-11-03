#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 7
#SBATCH -t 10-00:00:00
#SBATCH -J gemma-wrapper
#SBATCH -e gemma-wrapper_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gemma-wrapper_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ

PATH=$PATH:/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools

ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/gemma_wrapper
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/autosomes_dap_phased_filt.vcf.gz

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#############################################################

mkdir -p $OUTNAME

SUBSET=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_samples.txt
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno.txt
COVARIATES=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_weight_covariate.txt

LD_PRUNED_GENO=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/${OUTNAME}/${CHR}_${OUTNAME}_pruned_noNA.geno

gemma-wrapper --cache-dir ./${OUTNAME}/${OUTNAME}_gemma_cache --json -- \
    -g $LD_PRUNED_GENO \
    -p $PHENO \
    -gk \
    -debug > ${OUTNAME}/${OUTNAME}_K.json

##In theory everything above this only needs to be run once per subset of samples

##Select single phenotype to run. 
NUM_PHENO=19
cut -f $NUM_PHENO $PHENO > ${OUTNAME}/${OUTNAME}_${NUM_PHENO}.txt

GENO=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/${OUTNAME}/${CHR}_${OUTNAME}.geno
MAP=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/autosomes_merge_map.txt
gemma-wrapper --input ${OUTNAME}/${OUTNAME}_K.json --permutate 50 --permute-phenotype ${OUTNAME}/${OUTNAME}_${NUM_PHENO}.txt -- \
  -g $GENO \
  -c $COVARIATES \
  -a $MAP \
  -lmm 2 \
  -debug > ${OUTNAME}/${OUTNAME}_${NUM_PHENO}_GWA_50.txt