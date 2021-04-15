#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J bcftools_call
#SBATCH -e bcftools_call_%A_%a.err
#SBATCH -o bcftools_call_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools
cd $WORK_D

FILELIST=`ls -1 splits/split_* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
#was 200 per line, now 100
#25 lines .submit 24 bc 25th ran

##this could be run without 01, as it just checks if the file exists.
##the node crashed for some files

while read FILENAME
do
  FILE=$(basename $FILENAME)
  SAMPLE=${FILE/.sort.bam/}

  BAM=$FILENAME
  VCF=reference_panel/finch_324samps.postrecomb.sites.vcf.gz
  TSV=reference_panel/finch_324samps.postrecomb.sites.tsv.gz
  REFGEN=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
  OUT=${SAMPLE}.vcf.gz
if [ -f ${SAMPLE}.vcf.gz ]
  then
    echo $SAMPLE "exists"
  else
    #echo $SAMPLE "doesnt exist"
    bcftools mpileup --threads 2 -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
    bcftools index -f ${OUT}
  fi
done < $FILELIST
