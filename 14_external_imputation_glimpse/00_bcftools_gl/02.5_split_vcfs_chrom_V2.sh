#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 4-00:00:00
#SBATCH -J merge_bcftools
#SBATCH -e merge_bcftools_%A_%a.err
#SBATCH -o merge_bcftools_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools

#new 2020 data
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_round2_data
cd $WORK_D

# ##do once
# mkdir -p splits_vcfs
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

#when rerunning:
# realpath *vcf.gz > all_vcf_files.txt
# split -l 100 -d all_vcf_files.txt splits_vcfs/split_

FILELIST=`ls -1 splits_vcfs/split_* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
INTERVALS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/list_of_all_chromosomes.txt
#12 splits

mkdir -p merged_intervals

while read INTERVAL
do
  NAME=$(basename $FILELIST)
  bcftools merge -r $INTERVAL -l $FILELIST -O z -o merged_intervals/${NAME}_${INTERVAL}.vcf.gz
  bcftools index merged_intervals/${NAME}_${INTERVAL}.vcf.gz
done < $INTERVALS
