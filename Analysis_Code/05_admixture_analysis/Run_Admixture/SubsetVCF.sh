#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J SubsetVCF
#SBATCH -M snowy

module load bioinfo-tools vcftools
vcftools --gzvcf /proj/snic2020-2-19/private/darwins_finches/users/ashsende/Admixture_Oct21/resources/autosomes_all_highcov_lowcov_phased.vcf.gz \
--keep projAdmix.list --maf 0.05 --max-missing 1 \
--out autosomes_all_highcov_lowcov_phased.Daphne4spp.minMAF0.01.noMissing --recode
