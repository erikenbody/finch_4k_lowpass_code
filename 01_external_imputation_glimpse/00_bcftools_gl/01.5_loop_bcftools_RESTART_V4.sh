#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5-00:00:00
#SBATCH -J bcftools_call
#SBATCH -e bcftools_call_%A_%a.err
#SBATCH -o bcftools_call_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_mered_data
cd $WORK_D

#these are the duplicated samples that need to be re run as the merged bam file
realpath /crex/proj/snic2020-2-19/private/darwins_finches/alignment/lowpass/Fall2020_TH-2631/merged_duplicate_samples/*sorted.bam > duplicated_bams.txt

FILELIST=duplicated_bams.txt
REF_PANEL=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/reference_panel

while read FILENAME
do
  FILE=$(basename $FILENAME)
  SAMPLE=${FILE/.sorted.bam/}

  BAM=$FILENAME
  VCF=$REF_PANEL/finch_324samps.postrecomb.sites.vcf.gz
  TSV=$REF_PANEL/finch_324samps.postrecomb.sites.tsv.gz
  REFGEN=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
  OUT=${SAMPLE}.vcf.gz
if [ -f ${SAMPLE}.vcf.gz ]
  then
    echo $SAMPLE "exists"
  else
    #echo $SAMPLE "doesnt exist"
    bcftools mpileup --threads 4 -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
    bcftools index -f ${OUT}
  fi
done < $FILELIST
