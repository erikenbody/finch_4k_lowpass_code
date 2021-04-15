#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 10
#SBATCH -M rackham
#SBATCH -t 0-36:00:00
#SBATCH -J glimpse
#SBATCH -e glimpse_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o glimpse_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10

#V1.10 install Dec 7
GLIMPSE_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/GLIMPSE_1.10/GLIMPSE
#WORK_D=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/lowpass/external_imputation_glimpse
#cd $WORK_D

#with shapeit4 phased ref
SUBSET=${1}
REF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/${SUBSET}.324samps.postrecomb.bcf
MAP_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/13.5_shapeit4_phasing/postrecomb/gmap_files

#25 intervals for original dataset
#INTERVAL_FILES=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/merged_intervals
#13 intervals for 2020 lowpass data. Now 12 in march
#INTERVAL_FILES=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_round2_data/merged_intervals
#6 intervals feb 2021
INTERVAL_FILES=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/lowpass_round3_data/merged_intervals

VCF=`ls -1 $INTERVAL_FILES/*${SUBSET}.vcf.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
##FOR TESTING###
##VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/calc_gl_bcftools/merged_intervals/split_24_vcf.list_chr3.vcf.gz
####
#indexing now in other path
##bcftools index -f $VCF
mkdir -p $SUBSET
cd $SUBSET
#CHUNK. this changed vs. v 1.10 to only take the reference panel
$GLIMPSE_DIR/chunk/bin/GLIMPSE_chunk --input $REF --window-size 2000000 --buffer-size 200000 --region ${SUBSET} --output chunks.${SUBSET}_${SLURM_ARRAY_TASK_ID}.txt

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
