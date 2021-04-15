#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10-00:00:00
#SBATCH -J bcftools_call
#SBATCH -e bcftools_call_%A_%a.err
#SBATCH -o bcftools_call_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_round3_data
cd $WORK_D

#erik run and comment out when ready
# realpath /crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/map_new_reads_sentieon_2021/*/*sorted.bam > all_bam_files.txt
# split -l 100 -d all_bam_files.txt split_
# mkdir -p splits
# mv split_* splits
## should be six files


FILELIST=`ls -1 splits/split_* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
REF_PANEL=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/reference_panel

# while read FILENAME
# do
#   FILE=$(basename $FILENAME)
#   SAMPLE=${FILE/.sorted.bam/}
#
#   BAM=$FILENAME
#   VCF=$REF_PANEL/finch_324samps.postrecomb.sites.vcf.gz
#   TSV=$REF_PANEL/finch_324samps.postrecomb.sites.tsv.gz
#   REFGEN=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
#   OUT=${SAMPLE}.vcf.gz
# if [ -f ${SAMPLE}.vcf.gz ]
#   then
#     echo $SAMPLE "exists"
#   else
#     #echo $SAMPLE "doesnt exist"
#     bcftools mpileup --threads 4 -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
#     bcftools index -f ${OUT}
#   fi
# done < $FILELIST

#rerun one sample
FILENAME=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/map_new_reads_sentieon_2021/00Dap17863/00Dap17863.sorted.bam
FILE=$(basename $FILENAME)
SAMPLE=${FILE/.sorted.bam/}

BAM=$FILENAME
VCF=$REF_PANEL/finch_324samps.postrecomb.sites.vcf.gz
TSV=$REF_PANEL/finch_324samps.postrecomb.sites.tsv.gz
REFGEN=/crex/proj/snic2020-2-19/private/darwins_finches/reference/Camarhynchus_parvulus_V1.0.fasta
OUT=${SAMPLE}.vcf.gz
bcftools mpileup --threads 4 -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}
