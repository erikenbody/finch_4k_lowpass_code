#!/bin/bash
#SBATCH -A snic2020-15-109
#SBATCH -p core -n 5
#SBATCH -M snowy
#SBATCH -t 0-06:00:00
#SBATCH -J glimpse
#SBATCH -e glimpse_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o glimpse_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bcftools/1.10 samtools/1.10
ml python/2.7.15

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/test_sample
cd $WORK_D

# #one of the highest coverage samples
# BAM=/crex/proj/snic2020-2-19/private/darwins_finches/alignment/2017_Alignment/fortis_Daphne_5560.MD.RG.bam
# samtools view -s 1.05 -bo fortis_Daphne_5560_chr28.bam $BAM chr28
# samtools index $WORK_D/fortis_Daphne_5560_chr28.bam
# #check size of a high coverage bam
# samtools view -bo fortis_Daphne_5560_chr28_full $BAM chr28
# #check size of a normal low cov sample
# samtools view -bo 01Dap19399_chr28.bam /crex/proj/snic2020-2-19/private/darwins_finches/alignment/lowpass/01Dap19399.sort.bam chr28

#the subsampled one is a little bit higher than 19399 (6.7M vs 4.6M). Sounds about right.

sh ~/fc/finch_code/lowpass/14_external_imputation_glimpse/test_sample/lowpass_sentieon.sh fortis_Daphne_5560_chr28 $WORK_D/${CHR}/fortis_Daphne_5560_${CHR}.bam

#the glimpse script I wrote, but replaced the intervals with 5560 gvcf file
sbatch -p devel -t 00:10:00 -n 10 --array=1 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/test_sample/run_glimpse.sh chr28

#runs fast

VCF=GLIMPSE_ligated/chr28_1.merged.bcf
REF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/chr28.324samps.postrecomb.bcf
bcftools view -s fortis_Daphne_5560 $REF -O b -o fortis_Daphne_5560_reference.bcf


bcftools view -s fortis_Daphne_5560 finches_sentieon_27_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz -O b -o fortis_Daphne_5560_reference_not_phased.bcf

echo chr28 fortis_Daphne_5560_reference_not_phased.bcf fortis_Daphne_5560_reference_not_phased.bcf GLIMPSE_ligated/chr28_1.merged.bcf > concordance.lst
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE
module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
mkdir -p GLIMPSE_concordance

$GLIMPSE_DIR/concordance/bin/GLIMPSE_concordance --input concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 \
--bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000

ml matplotlib/3.0.3-foss-2019a-Python-3.7.2
#for running their plotting script
cd GLIMPSE_concordance/output
/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE/tutorial/plot/concordance_plot.py
