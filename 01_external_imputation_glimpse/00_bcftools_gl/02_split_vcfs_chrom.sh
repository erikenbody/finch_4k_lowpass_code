#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 4-00:00:00
#SBATCH -J split_hc
#SBATCH -e split_hc_%A_%a.err
#SBATCH -o split_hc_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools
cd $WORK_D

#I ran this for the 25th array already as a test

# ##do once
mkdir -p splits_vcfs
# for FILELIST in splits/split_*
# do
#   NAME2=$(basename $FILELIST)
#   while read FILENAME
#     do
#       FILE=$(basename $FILENAME)
#       SAMPLE=${FILE/.sort.bam/}
#       realpath ${SAMPLE}.vcf.gz >> splits_vcfs/${NAME2}_vcf.list
#     done < $FILELIST
# done

FILELIST=`ls -1 splits_vcfs/split_* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
INTERVALS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/chr_list.txt

while read FILE
do
  #ls ${FILE}.csi
  bcftools index $FILE
done < $FILELIST

mkdir -p merged_intervals

while read INTERVAL
do
  NAME=$(basename $FILELIST)
  bcftools merge -r $INTERVAL -l $FILELIST -O z -o merged_intervals/${NAME}_${INTERVAL}.vcf.gz
  bcftools index merged_intervals/${NAME}_${INTERVAL}.vcf.gz
done < $INTERVALS
