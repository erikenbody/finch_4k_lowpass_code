#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -M rackham
#SBATCH -p devel
#SBATCH -n 1
#SBATCH -t 0-00:10:00
#SBATCH -J split_hc
#SBATCH -e split_hc_%A_%a.err
#SBATCH -o split_hc_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_round2_data
cd $WORK_D

#I ran this for the 25th array already as a test

# ##do once
mkdir -p splits_vcfs
for FILELIST in splits/split_*
do
  NAME2=$(basename $FILELIST)
  while read FILENAME
    do
      FILE=$(basename $FILENAME)
      SAMPLE=${FILE/.sorted.bam/}
      realpath ${SAMPLE}.vcf.gz >> splits_vcfs/${NAME2}_vcf.list
    done < $FILELIST
done
