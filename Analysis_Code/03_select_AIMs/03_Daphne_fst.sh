#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 10
#SBATCH -t 2-00:00:00
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

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/pixy
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/lowcoverage_highcoverage_merge/autosomes_all_highcov_lowcov_phased.vcf.gz
SUBSET=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/24_HMGA2_ALX1/Daphne_geospiza_n50.txt
#grep -v "Bar" $SUBSET > just_Daphne_n50.txt
#bcftools view -S just_Daphne_n50.txt -O z -o Daphne_n50.vcf.gz $VCF
#tabix Daphne_n50.vcf.gz
POPS_FILE=Daphne_geospiza_n50_pops.txt

pixy --stats fst \
  --vcf Daphne_n50.vcf.gz \
  --window_size 1 \
  --n_cores 5 \
  --populations $POPS_FILE \
  --bypass_invariant_check yes \
  --output_folder output_Daphne \
  --output_prefix output_Daphne_persite
