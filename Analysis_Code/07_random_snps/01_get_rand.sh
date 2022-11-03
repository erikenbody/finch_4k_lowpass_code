#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 0-12:00:00
#SBATCH -J getrand
#SBATCH -e getrand_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o getrand_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.14 BEDTools/2.29.2
ml python/2.7.15 FastTree/2.1.8

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/35_rand_snps
cd $WORK_D
SUBSET=autosomes
VCF_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/lowcoverage_highcoverage_merge
GZVCF=${VCF_D}/${SUBSET}_all_highcov_lowcov_phased.vcf.gz


BED1=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny/autosomes_for_ful_mag_species_lmm_PEAKS.bed
BED2=autosomes_Daphne_cluster1_multivariate_INTs_lmm1_PEAKS.bed
cat $BED1 $BED2 > combined_peaks.bed
DROP_BED=combined_peaks.bed

bedtools intersect -v -a $GZVCF -b $DROP_BED -wa -header > nonpeak_all_highcov_lowcov_phased.vcf
bgzip nonpeak_all_highcov_lowcov_phased.vcf
bcftools index nonpeak_all_highcov_lowcov_phased.vcf.gz

bcftools query -f '%CHROM\t%POS\n' nonpeak_all_highcov_lowcov_phased.vcf.gz | wc -l > num_sites.txt
VAL=`cat num_sites.txt`
SITES=`echo $((VAL / 100000))`
bcftools +prune -w $SITES -n 1 -N rand -O z -o nonpeak_all_highcov_lowcov_phased_pruned.vcf.gz nonpeak_all_highcov_lowcov_phased.vcf.gz
