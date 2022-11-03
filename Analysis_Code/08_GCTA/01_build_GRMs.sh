#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 20
#SBATCH -t 1-00:00:00
#SBATCH -J gcta
#SBATCH -e gcta_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gcta_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools qctool/2.0.6 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/37_GCTA
cd $WORK_D

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
#OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
OUTNAME=Daphne_four_species

#############################################################

#rerun with make-grm-alg 1
gcta64 --mbfile ${OUTNAME}/path_to_autosomes_bfiles.txt --make-grm-alg 1 --maf 0.05 --make-grm-gz --thread-num 8 --out ${OUTNAME}/${CHR}_${OUTNAME}_auto_flat_maf05_algo1 --autosome-num 49
OUTNAME=Daphne_cluster1
gcta64 --mbfile ${OUTNAME}/path_to_autosomes_bfiles.txt --make-grm-alg 1 --maf 0.05 --make-grm-gz --thread-num 8 --out ${OUTNAME}/${CHR}_${OUTNAME}_auto_flat_maf05_algo1 --autosome-num 49
