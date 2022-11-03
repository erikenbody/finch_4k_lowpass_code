#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 5
#SBATCH -t 1-00:00:00
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
INTERVALS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/full_ordered_chr.list

mkdir -p $OUTNAME/per_chr
mkdir -p results_${OUTNAME}/per_chr
while read chr
do
  gcta64 --bfile ${OUTNAME}/${chr}_${OUTNAME} --make-grm --threads 5 --out ${OUTNAME}/per_chr/${chr}_${OUTNAME} --autosome-num 49
  gcta64 --reml --grm ${OUTNAME}/per_chr/${chr}_${OUTNAME} --pheno $PHENO --out results_${OUTNAME}/per_chr/${chr}_reml_PC1
done < $INTERVALS

cat results_Daphne_cluster1/per_chr/*hsq | grep "V(G)/Vp" | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'
cat results_Daphne_cluster1/per_chr/*hsq | grep "V(G)/Vp" | awk '{ sum += $2; n++ } END { if (n > 0) print sum ;}'

while read chr
do
  echo ${OUTNAME}/per_chr/${chr}_${OUTNAME} >> $OUTNAME/per_chr_grms.txt
done < $INTERVALS

grep -v "chrZ" $OUTNAME/per_chr_grms.txt > $OUTNAME/autosomes_per_chr_grms.txt

gcta64 --reml --mgrm $OUTNAME/autosomes_per_chr_grms.txt --pheno $PHENO --out results_${OUTNAME}/compare_all_chr_greml

grep /Vp results_${OUTNAME}/compare_all_chr_greml.hsq | grep -v "Sum" > results_${OUTNAME}/compare_all_chr_greml_VGVP.txt

gcta64 --reml --mgrm $OUTNAME/per_chr_grms.txt --pheno $PHENO --out results_${OUTNAME}/compare_all_chr_greml_incluZ

grep /Vp results_${OUTNAME}/compare_all_chr_greml_incluZ.hsq | grep -v "Sum" > results_${OUTNAME}/compare_all_chr_greml_VGVP_incluZ.txt
