#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 8
#SBATCH -t 2-00:00:00
#SBATCH -J grml
#SBATCH -e grml_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o grml_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com


ml bioinfo-tools qctool/2.0.6 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/37_GCTA
cd $WORK_D

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#############################################################

PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_gcta.txt
COV=~/fc/finch_code/lowpass/37_GCTA/${OUTNAME}_qtl_covs.txt

gcta64 --reml --grm ${OUTNAME}/snp_group1 --covar $COV --pheno $PHENO --out results_${OUTNAME}/snp_group1_greml_with_qtls.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --pheno $PHENO --out results_${OUTNAME}/greml_PC1_maf05.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --covar $COV --pheno $PHENO --out results_${OUTNAME}/greml_PC1_with_qtls_maf05.txt

#just weight
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_weight_INTs.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --covar $COV --pheno $PHENO --out results_${OUTNAME}/greml_weight_INts_with_qtls.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --pheno $PHENO --out results_${OUTNAME}/greml_weight_INts.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --pheno $PHENO --out results_${OUTNAME}/greml_weight_INts_maf05.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --covar $COV --pheno $PHENO --out results_${OUTNAME}/greml_weight_INts_with_qtls_maf05.txt
#PC2
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_PC2_gcta.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto --covar $COV --pheno $PHENO --out results_${OUTNAME}/greml_PC2_with_qtls.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --covar $COV --pheno $PHENO --out results_${OUTNAME}/greml_PC2_with_qtls_maf05.txt
gcta64 --reml --grm ${OUTNAME}/${CHR}_${OUTNAME}_auto_maf05 --pheno $PHENO --out results_${OUTNAME}/greml_PC2_maf05.txt
