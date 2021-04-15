#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 10
#SBATCH -M rackham
#SBATCH -t 0-36:00:00
#SBATCH -J glimpse
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10

#V1.10 install Dec 7
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE
#WORK_D=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/lowpass/external_imputation_glimpse
#cd $WORK_D

#with shapeit4 phased ref
SUBSET=${1}
REF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/${SUBSET}.324samps.postrecomb.bcf
bcftools view -s ^fortis_Daphne_5560 -O b -o ${SUBSET}.323samps.postrecomb.bcf $REF
REF=${SUBSET}.323samps.postrecomb.bcf
bctools index $REF

MAP_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/gmap_files

INTERVAL_FILES=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/sentieon_haplotyper/merged_intervals
#VCF=`ls -1 $INTERVAL_FILES/*${SUBSET}.vcf.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

#using subsampled high cov individual
#VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/test_sample/fortis_Daphne_5560_chr28/fortis_Daphne_5560_chr28.g.vcf.gz
VCF=${2}
bcftools index -f $VCF

#CHUNK. this changed vs. v 1.10 to only take the reference panel
$GLIMPSE_DIR/chunk/bin/GLIMPSE_chunk --reference $REF --input $VCF --window-size 2000000 --buffer-size 200000 --region ${SUBSET} --output chunks.${SUBSET}_${SLURM_ARRAY_TASK_ID}.txt

mkdir -p GLIMPSE_imputed

#PHASE, note that 16 iterations appears to be the max number. After 1.1 15 is max
while IFS="" read -r LINE || [ -n "$LINE" ];
do
  printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
  IRG=$(echo $LINE | cut -d" " -f3)
  ORG=$(echo $LINE | cut -d" " -f4)
  OUT=GLIMPSE_imputed/merged.${SUBSET}_${SLURM_ARRAY_TASK_ID}.${ID}.bcf
  $GLIMPSE_DIR/phase/bin/GLIMPSE_phase --main 15 --map $MAP_DIR/${SUBSET}_parv.gmap --input ${VCF} --reference ${REF} --input-region ${IRG} --output-region ${ORG} --output ${OUT} --thread 10
  bcftools index -f ${OUT}
done < chunks.${SUBSET}_${SLURM_ARRAY_TASK_ID}.txt

mkdir -p GLIMPSE_ligated

#LIGATE
LST=GLIMPSE_ligated/${SUBSET}_${SLURM_ARRAY_TASK_ID}.txt
ls GLIMPSE_imputed/merged.${SUBSET}_${SLURM_ARRAY_TASK_ID}.*.bcf > ${LST}
OUT=GLIMPSE_ligated/${SUBSET}_${SLURM_ARRAY_TASK_ID}.merged.bcf
$GLIMPSE_DIR/ligate/bin/GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}

mkdir -p GLIMPSE_phased

#SAMPLE
VCF=GLIMPSE_ligated/${SUBSET}_${SLURM_ARRAY_TASK_ID}.merged.bcf
OUT=GLIMPSE_phased/${SUBSET}_${SLURM_ARRAY_TASK_ID}.phased.bcf
$GLIMPSE_DIR/sample/bin/GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT}

#format
# cd GLIMPSE_phased
# VCF=${SUBSET}_${SLURM_ARRAY_TASK_ID}.phased.bcf
# SAMPLE=${SUBSET}_${SLURM_ARRAY_TASK_ID}
# bcftools query -l $VCF | tr "\n" "\t" > ${SAMPLE}_header.txt
# echo "" >> ${SAMPLE}_header.txt #add new line at end of header
# bcftools query -f '[%GT\t]\n' $VCF > ${SAMPLE}_raw.genos
# cat ${SAMPLE}_header.txt ${SAMPLE}_raw.genos > ${SAMPLE}.genos
# bcftools query -f '%CHROM\t %POS\t %REF\t %ALT\n' $VCF > ${SAMPLE}_raw.pos
# echo -e "CHROM\tPOS\tREF\tALT" | cat - ${SAMPLE}_raw.pos > ${SAMPLE}.pos
# paste -d "\t" ${SAMPLE}.pos ${SAMPLE}.genos > ${SAMPLE}.pos.genos
# rm *_raw*

# bcftools view -r chr1A -v snps -O z -o ${SAMPLE}_SNPs_raw.genos.vcf.gz $VCF
# VCF=01Dap18078_SNPs_raw.genos.vcf.gz
# bcftools query -l $VCF | tr "\n" "\t" > ${SAMPLE}_header.txt
# echo "" >> ${SAMPLE}_header.txt #add new line at end of header
# bcftools query -f '[%GT\t]\n' $VCF > ${SAMPLE}_raw.genos
# cat ${SAMPLE}_header.txt ${SAMPLE}_raw.genos > ${SAMPLE}.genos
# bcftools query -f '%CHROM\t %POS\t %REF\t %ALT\n' $VCF > ${SAMPLE}_raw.pos
# echo -e "CHROM\tPOS\tREF\tALT" | cat - ${SAMPLE}_raw.pos > ${SAMPLE}.pos
# paste -d "\t" ${SAMPLE}.pos ${SAMPLE}.genos > ${SAMPLE}.pos.genos
# rm *_raw*
