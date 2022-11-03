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

##run with --make-grm-alg 1 which should better handle off diagonals that have GRM > 1 due to rare alleles

for group in snp_group_median1 snp_group_median2
do
  gcta64 --make-grm-alg 1 --maf 0.05 --threads 3 --mbfile ${OUTNAME}/path_to_autosomes_bfiles.txt --extract ${OUTNAME}/${OUTNAME}_${group}.txt --make-grm --out ${OUTNAME}/${group}_maf05_2ldgroups_V2 --autosome-num 49
  echo ${OUTNAME}/${group}_maf05_2ldgroups_V2 >> ${OUTNAME}/multiple_GRMs_ld_maf05_2ldgroups_V2.txt
done

#set a max bin

for group in snp_group_median1 snp_group_median2
do
  gcta64 --make-grm-alg 1 --max-maf 0.05 --threads 3 --mbfile ${OUTNAME}/path_to_autosomes_bfiles.txt --extract ${OUTNAME}/${OUTNAME}_${group}.txt --make-grm --out ${OUTNAME}/${group}_maxmaf05_2ldgroups_V2 --autosome-num 49
  echo ${OUTNAME}/${group}_maxmaf05_2ldgroups_V2 >> ${OUTNAME}/multiple_GRMs_ld_maxmaf05_2ldgroups_V2.txt
done

cat ${OUTNAME}/multiple_GRMs_ld_maf05_2ldgroups_V2.txt ${OUTNAME}/multiple_GRMs_ld_maxmaf05_2ldgroups_V2.txt > ${OUTNAME}/multiple_GRMs_two_maf_bins_two_median_ld_bins_V2.txt

gcta64 --reml --threads 3 --mgrm ${OUTNAME}/multiple_GRMs_two_maf_bins_two_median_ld_bins_V2.txt --pheno $PHENO --out results_${OUTNAME}/ldms_reml_results_2mafbins_2median_ldbins_V2.txt

#run with 20 PCs
gcta64 --reml --threads 3 --mgrm ${OUTNAME}/multiple_GRMs_two_maf_bins_two_median_ld_bins_V2.txt --qcovar $OUTNAME/PCA_autosomes_Daphne_cluster1_auto_maf05.eigenvec --pheno $PHENO --out results_${OUTNAME}/ldms_reml_results_2mafbins_2median_ldbins_20PCcov_V2.txt


#run for weight and PC2
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_PC2_gcta.txt
gcta64 --reml --threads 3 --mgrm ${OUTNAME}/multiple_GRMs_two_maf_bins_two_median_ld_bins_V2.txt --pheno $PHENO --out results_${OUTNAME}/bill_PC2_ldms_reml_results_2mafbins_2median_ldbins_V2.txt

PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno_weightINT_gcta.txt
gcta64 --reml --threads 3 --mgrm ${OUTNAME}/multiple_GRMs_two_maf_bins_two_median_ld_bins_V2.txt --pheno $PHENO --out results_${OUTNAME}/weightINT_ldms_reml_results_2mafbins_2median_ldbins_V2.txt
