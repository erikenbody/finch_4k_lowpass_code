#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4 #max mem when not using pruned dataset
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J related
#SBATCH -e related_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o related_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ
ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/autosomes_dap_phased_filt_pruned_20kb.vcf.gz
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/sample_pheno_all_3957_pheno.txt

CHR=autosomes

##formatting of map and SNP id only ever needs to be run once
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o ${CHR}_merge_ID_pruned.vcf.gz $VCF
##this makes a bit of a mess, but basically the problem is that qctools which makes the bimbim genotype chr_format
##outputs a snp ID as ID:chrom:pos. And the only ID thats easy to make is chrom_ID_pos_ref_alt
##could make it prettier by scripting my own qctool like thing, but why bother?
bcftools query -f '%ID\:%CHROM\:%POS, %POS\, %CHROM\n' ${CHR}_merge_ID_pruned.vcf.gz > ${CHR}_merge_map_pruned.txt

qctool -g ${CHR}_merge_ID_pruned.vcf.gz -ofiletype bimbam_dosage -og ${CHR}.geno

##create relatedness matrix (centered is best reccomendation for starting off in the documentation)

grep -v "NA" ${CHR}.geno  > ${CHR}_noNA.geno

$GEMMA/gemma -g ${CHR}_noNA.geno -p $PHENO -gk 1 -o ${CHR}_noNA.relate

mkdir -p Daphne_3957
cp output/* Daphne_3957
