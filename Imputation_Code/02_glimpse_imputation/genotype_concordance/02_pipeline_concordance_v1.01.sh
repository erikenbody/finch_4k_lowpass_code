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
#using old install of glimpse
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/test_sample
cd $WORK_D

#23654

CHR=chr12
mkdir -p $CHR
cd $CHR

#one of the highest coverage samples
BAM=/crex/proj/snic2020-2-19/private/darwins_finches/alignment/2017_Alignment/fortis_Daphne_5560.MD.RG.bam
samtools view -s 1.05 -bo fortis_Daphne_5560_${CHR}.bam $BAM $CHR
samtools index fortis_Daphne_5560_${CHR}.bam
ml python/2.7.15
sh ~/fc/finch_code/lowpass/14_external_imputation_glimpse/test_sample/lowpass_sentieon.sh fortis_Daphne_5560_${CHR} $WORK_D/${CHR}/fortis_Daphne_5560_${CHR}.bam

GVCF_PATH=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/test_sample/fortis_Daphne_5560_${CHR}/fortis_Daphne_5560_${CHR}.g.vcf.gz

#the glimpse script I wrote, but replaced the intervals with 5560 gvcf file

sbatch -p devel -t 00:10:00 -n 10 --array=1 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/test_sample/run_glimpse_V1.01.sh $CHR ${GVCF_PATH}

#runs fast

VCF=GLIMPSE_ligated/${CHR}_1.merged.bcf
#REF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/${CHR}.324samps.postrecomb.bcf
#bcftools view -s fortis_Daphne_5560 $REF -O b -o fortis_Daphne_5560_reference.bcf

OGVCF=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/07_BIALLELIC_ONLY/finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz

bcftools view -s fortis_Daphne_5560 $OGVCF -O b -o fortis_Daphne_5560_reference_not_phased_${CHR}.bcf $CHR

echo ${CHR} fortis_Daphne_5560_reference_not_phased_${CHR}.bcf fortis_Daphne_5560_reference_not_phased_${CHR}.bcf GLIMPSE_ligated/${CHR}_1.merged.bcf > concordance.lst
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE
module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
mkdir -p GLIMPSE_concordance

$GLIMPSE_DIR/concordance/bin/GLIMPSE_concordance --input concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 \
--bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000

ml matplotlib/3.0.3-foss-2019a-Python-3.7.2
#for running their plotting script
cd GLIMPSE_concordance/
/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE/tutorial/plot/concordance_plot.py
