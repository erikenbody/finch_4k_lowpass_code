#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 16
#SBATCH -t 1-00:00:00
#SBATCH -J pixy
#SBATCH -e pixy_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o pixy_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml conda/latest bcftools/1.12
source conda_init.sh
CONDA_ENVS_PATH=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools
#conda install -c conda-forge pixy
#conda create --name pixy -c conda-forge pixy
#to updatee
#conda install --yes -c conda-forge pixy
conda activate pixy
INTERVAL_LIST=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPLES=four_species_for_AIM_samples.txt

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/pixy
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/lowcoverage_highcoverage_merge/autosomes_all_highcov_lowcov_phased.vcf.gz
#bcftools view -r $INTERVAL -S $SAMPLES -O z -o four_species_${INTERVAL}.vcf.gz $VCF
tabix four_species_${INTERVAL}.vcf.gz
POPS_FILE=four_species_for_AIM.txt

pixy --stats fst \
  --vcf four_species_${INTERVAL}.vcf.gz \
  --chromosomes $INTERVAL \
  --window_size 1 \
  --n_cores 8 \
  --populations $POPS_FILE \
  --bypass_invariant_check yes \
  --output_folder output_PERSITE \
  --output_prefix four_species_pixy_${INTERVAL}_persite
