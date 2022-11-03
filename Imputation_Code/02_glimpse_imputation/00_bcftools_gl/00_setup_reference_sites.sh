#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH -J ref_sites
#SBATCH -e ref_sites_%A_%a.err
#SBATCH -o ref_sites_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools
cd $WORK_D
OGVCF=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/07_BIALLELIC_ONLY/finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz

mkdir -p reference_panel
bcftools view -G -m 2 -M 2 -v snps $OGVCF -Oz -o reference_panel/finch_324samps.postrecomb.sites.vcf.gz
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/finch_324samps.postrecomb.sites.vcf.gz | bgzip -c > reference_panel/finch_324samps.postrecomb.sites.tsv.gz
tabix -s1 -b2 -e2 reference_panel/finch_324samps.postrecomb.sites.tsv.gz
