#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 10
#SBATCH -M rackham
#SBATCH -t 00-9:00:00 #may need to up this
#SBATCH -J snmf
#SBATCH -e snmf_%A_%a.err            # File to which STDERR will be written
#SBATCH -o snmf_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools plink2/2.00-alpha-2.3-20200124 bcftools/1.12 R_packages/4.0.0

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/snmf/aims
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/lowcoverage_highcoverage_merge/autosomes_all_highcov_lowcov_phased.vcf.gz

##Using strict neutral filtering  (i.e. no chromosomes with peak from Sci Adv paper)

AIMS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/pixy/autosomal_four_species_V3/autosomal_four_species_V3_aims_neutral_STRICT.txt
OUTNAME=Daphne_neutral_aims_strict

#bcftools view -R $AIMS -O z -o ${OUTNAME}_Daphne_aims.vcf.gz $VCF
#prune for linked sites
bcftools +prune -m 0.3 -w 20kb --nsites-per-win-mode rand --random-seed 42 ${OUTNAME}_Daphne_aims.vcf.gz -O v -o ${OUTNAME}_pruned_Daphne_AIMs.vcf
#get sample names
bcftools query -l ${OUTNAME}_pruned_Daphne_AIMs.vcf > ${OUTNAME}_pruned_Daphne_AIMs_samples.txt

Rscript ~/fc/finch_code/lowpass/22_snmf/run_snmf_cluster_Daphne.R ${OUTNAME}_pruned_Daphne_AIMs.vcf ${OUTNAME}_pruned_Daphne_AIMs_samples.txt $OUTNAME 10

##using less strict neutral filtering (i.e. only peaks removed)

AIMS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/pixy/autosomal_four_species_V3/autosomal_four_species_V3_aims_neutral.txt
OUTNAME=Daphne_neutral_aims

b#cftools view -R $AIMS -O z -o ${OUTNAME}_Daphne_aims.vcf.gz $VCF
#prune for linked sites
bcftools +prune -m 0.3 -w 20kb --nsites-per-win-mode rand --random-seed 42 ${OUTNAME}_Daphne_aims.vcf.gz -O v -o ${OUTNAME}_pruned_Daphne_AIMs.vcf
#get sample names
bcftools query -l ${OUTNAME}_pruned_Daphne_AIMs.vcf > ${OUTNAME}_pruned_Daphne_AIMs_samples.txt

Rscript ~/fc/finch_code/lowpass/22_snmf/run_snmf_cluster_Daphne.R ${OUTNAME}_pruned_Daphne_AIMs.vcf ${OUTNAME}_pruned_Daphne_AIMs_samples.txt $OUTNAME 10
