#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 2-00:00:00
#SBATCH -J merge_glimpse
#SBATCH -e merge_glimpse_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o glimpse_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse
cd $WORK_D

#SUBSET=`ls -d */ | cut -f1 -d'/' | grep 'chr' | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
INTERVAL_LIST=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/full_ordered_chr.list
SUBSET=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

mkdir -p merged_vcfs
LST=${SUBSET}_ligated.list
ls ${SUBSET}/${SUBSET}/GLIMPSE_ligated/${SUBSET}*.bcf > merged_vcfs/${LST}
ls lowpass_round2_data/${SUBSET}/${SUBSET}/GLIMPSE_ligated/${SUBSET}*.bcf >> merged_vcfs/${LST}
ls lowpass_round3_data/${SUBSET}/${SUBSET}/GLIMPSE_ligated/${SUBSET}*.bcf >> merged_vcfs/${LST}

bcftools merge -l merged_vcfs/${LST} -O z -o merged_vcfs/${SUBSET}_ligated_merge.vcf.gz
bcftools index merged_vcfs/${SUBSET}_ligated_merge.vcf.gz


LST=${SUBSET}.list
ls ${SUBSET}/${SUBSET}/GLIMPSE_phased/${SUBSET}*.bcf > merged_vcfs/${LST}
ls lowpass_round2_data/${SUBSET}/${SUBSET}/GLIMPSE_phased/${SUBSET}*.bcf >> merged_vcfs/${LST}
ls lowpass_round3_data/${SUBSET}/${SUBSET}/GLIMPSE_phased/${SUBSET}*.bcf >> merged_vcfs/${LST}

bcftools merge -l merged_vcfs/${LST} -O z -o merged_vcfs/${SUBSET}_phased_merge.vcf.gz
bcftools index merged_vcfs/${SUBSET}_phased_merge.vcf.gz
