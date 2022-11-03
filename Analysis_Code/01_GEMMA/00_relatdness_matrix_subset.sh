#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 1 #dont need much to run this on the pruned dataset. interactive is fine
#SBATCH -t 1-00:00:00
#SBATCH -J related
#SBATCH -e related_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o related_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ

#ml bioinfo-tools GEMMA/0.98.1 vcftools/0.1.16  plink2/2.00-alpha-2-20190429 qctool/2.0.6 bcftools
ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/autosomes_dap_phased_filt_pruned_20kb.vcf.gz

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#OUTNAME=fortis_bill_color
#OUTNAME=scandens_bill_color

#############################################################

SUBSET=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_samples.txt
PHENO=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_pheno.txt

mkdir -p $OUTNAME

##Can start here if Daphne birds
qctool -g ${CHR}_merge_ID_pruned.vcf.gz -ofiletype bimbam_dosage -og $OUTNAME/${CHR}_${OUTNAME}_pruned.geno -incl-samples $SUBSET

##create relatedness matrix (centered is best reccomendation for starting off in the documentation)

grep -v "NA" $OUTNAME/${CHR}_${OUTNAME}_pruned.geno  > $OUTNAME/${CHR}_${OUTNAME}_pruned_noNA.geno

$GEMMA/gemma -g $OUTNAME/${CHR}_${OUTNAME}_pruned_noNA.geno -p $PHENO -gk 1 -o ${CHR}_${OUTNAME}.noNA.relate
##cant control location of output so do it like this
cp output/${CHR}_${OUTNAME}.noNA.relate* $OUTNAME
