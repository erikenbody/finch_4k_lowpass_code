#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 7
#SBATCH -M snowy
#SBATCH -t 2-00:00:00
#SBATCH -J whatshap
#SBATCH -e whatshap_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o whatshap_%J_%A_%a.out

module load bioinfo-tools bcftools/1.10 samtools/1.10
#whatshap installed with pip
export PATH=$HOME/.local/bin:$PATH

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/23_shapeit4_phasing_new_samps
cd $WORK_D

#31 chromosomes
INTERVAL_LIST=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt

#if you want to run it as a job array rather than a loop:
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
REF=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/07_BIALLELIC_ONLY/AC_FILTERED/finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC_ACf.vcf.gz

####
SAMPLE=${1}

##because I had renamed the VCF, need to track down the BAM file with the original name

MATCH_FILE=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/reheader_vcf_23Mar21.txt

##this is a disastor of organization, but it works! Sorry if you're trying to translate this in the future. Here is how I checked my work:
# while read SAMPLE
# do
#   OG_NAME=$(grep -hw $SAMPLE $MATCH_FILE | awk '{ print $1 }' FS='\t')
#   BAM_FILE=$(grep -hw ${OG_NAME}.deduped.bam all_bam_paths_messy.txt)
#   echo $BAM_FILE $SAMPLE
#   ls -1 $BAM_FILE >> test_bam_locs.txt
# done < vcf_sample_list.txt

OG_NAME=$(grep -hw $SAMPLE $MATCH_FILE | awk '{ print $1 }' FS='\t')
BAM_FILE=$(grep -hw ${OG_NAME}.deduped.bam all_bam_paths_messy.txt)

mkdir -p ${SNIC_TMP}/${SAMPLE}

bcftools view -r $INTERVAL -s $SAMPLE -O z -o ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz $VCF
bcftools index ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz

cd ${SNIC_TMP}/${SAMPLE}
scp $REF ${SNIC_TMP}/${SAMPLE}
scp $REF.fai ${SNIC_TMP}/${SAMPLE}

samtools view -b $BAM_FILE $INTERVAL > ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.MD.RG.bam
samtools index ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.MD.RG.bam

#ignore read groups is important here, because the VCF has been reheadered previously
whatshap phase --ignore-read-groups -o ${SAMPLE}_${INTERVAL}.whatshap.vcf --reference=Camarhynchus_parvulus_V1.0.fasta ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz ${SAMPLE}_${INTERVAL}.MD.RG.bam

mkdir -p ${WORK_D}/${SAMPLE}
mv ${SAMPLE}_${INTERVAL}.whatshap.vcf ${WORK_D}/${SAMPLE}

bgzip ${WORK_D}/${SAMPLE}/${SAMPLE}_${INTERVAL}.whatshap.vcf
bcftools index ${WORK_D}/${SAMPLE}/${SAMPLE}_${INTERVAL}.whatshap.vcf.gz
