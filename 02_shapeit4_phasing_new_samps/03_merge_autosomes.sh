#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 0-03:00:00
#SBATCH -J concat_glimpse
#SBATCH -e concat_glimpse_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o concat_glimpse_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bcftools/1.12

SUBSET=autosomes

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/23_shapeit4_phasing_new_samps
cd $WORK_D

LST=${SUBSET}.list
ls *postrecomb.bcf | grep -v "chrZ" > ${LST}

bcftools concat --threads 1 -f ${LST} -O z -o shapeit_${SUBSET}.vcf.gz
bcftools index shapeit_${SUBSET}.vcf.gz
