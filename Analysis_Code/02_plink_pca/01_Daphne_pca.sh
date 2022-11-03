#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 8
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J pca
#SBATCH -e pca_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o pca_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools vcftools zlib/1.2.11 plink2

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/17_plink_pca.sh
cd $WORK_D

#feb 2021, with new data
INPUT=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/autosomes_dap_phased_filt.vcf.gz
plink2 --vcf $INPUT --pca 4 --out Daphne_filtered_4ksamps_fixed --allow-extra-chr --const-fid --autosome-num 30 --threads 3
