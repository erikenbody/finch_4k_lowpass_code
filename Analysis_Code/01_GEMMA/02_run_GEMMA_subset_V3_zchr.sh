#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 3 #two needed after adding in new samps
#SBATCH -M rackham
#SBATCH -t 3-00:00:00
#SBATCH -J gemma
#SBATCH -e gemma_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gemma_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com
#SBATCH --cpus-per-task=1

ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/gemma-0.98.4
cd $WORK_D
CHR=chrZ

#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#############################################################

echo $OUTNAME

SUBSET=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_samples.txt
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno.txt
COVARIATES=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_weight_covariate.txt

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/chrZ_dap_phased_filt.vcf.gz

########this only really needs to be generated once. The vcf to be subsetted for input in gemma
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o $OUTNAME/${CHR}_merge_ID.vcf.gz $VCF
##this makes a bit of a mess, but basically the problem is that qctools which makes the bimbim genotype chr_format
##outputs a snp ID as ID:chrom:pos. And the only ID thats easy to make is chrom_ID_pos_ref_alt
##could make it prettier by scripting my own qctool like thing, but why bother?
bcftools query -f '%ID\:%CHROM\:%POS, %POS\, %CHROM\n' $OUTNAME/${CHR}_merge_ID.vcf.gz > $OUTNAME/${CHR}_merge_map.txt
###############################################################################################################
#this just gets run once per subset
qctool -g $OUTNAME/${CHR}_merge_ID.vcf.gz -ofiletype bimbam_dosage -og ${OUTNAME}/${CHR}_${OUTNAME}.geno -incl-samples $SUBSET
#######################################################################################################################

mkdir -p ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out

## INTs body weight with no beak size covariates, only sex and relatedness

$GEMMA/gemma -p $PHENO -c ~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_covariate.txt \
-a ${CHR}_merge_map.txt -g ${OUTNAME}/${CHR}_${OUTNAME}.geno \
  -k $OUTNAME/autosomes_${OUTNAME}.noNA.relate.cXX.txt -lmm 4 -outdir ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out \
  -n 31 \
  -o ${CHR}_${OUTNAME}_weight_INTs_no_beak_cov_lmm4

# ##or with the two PC.sp values
$GEMMA/gemma -p $PHENO -c  $COVARIATES \
  -a ${CHR}_merge_map.txt -g ${OUTNAME}/${CHR}_${OUTNAME}.geno \
  -k $OUTNAME/autosomes_${OUTNAME}.noNA.relate.cXX.txt -lmm 4 -outdir ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out \
  -n 24 25 \
  -o ${CHR}_${OUTNAME}_multivariate_PCsp_lmm4

##or with the two PC.sp values and no weight covariate
COVARIATES=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_covariate.txt
$GEMMA/gemma -p $PHENO -c  $COVARIATES \
  -a ${CHR}_merge_map.txt -g ${OUTNAME}/${CHR}_${OUTNAME}.geno \
  -k $OUTNAME/autosomes_${OUTNAME}.noNA.relate.cXX.txt -lmm 4 -outdir ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out \
  -n 24 25 \
  -o ${CHR}_${OUTNAME}_multivariate_PCsp_no_weight_lmm4

# run with top HMGA2 SNP as a covariate for PC.sp
COVARIATES=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_sex_weight_p33100090_covariate.txt
$GEMMA/gemma -p $PHENO -c  $COVARIATES \
  -a ${CHR}_merge_map.txt -g ${OUTNAME}/${CHR}_${OUTNAME}.geno \
  -k $OUTNAME/autosomes_${OUTNAME}.noNA.relate.cXX.txt -lmm 4 -outdir ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out \
  -n 24 \
  -o ${CHR}_${OUTNAME}_PC1sp_cov.p33100090_lmm4