#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J gcta
#SBATCH -e gcta_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gcta_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools qctool/2.0.6 bcftools/1.12

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/37_GCTA
cd $WORK_D

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
#OUTNAME=Daphne_cluster1
OUTNAME=Daphne_four_species
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
#############################################################

mkdir -p $OUTNAME

VCF_DIR=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs

while read chr
do
  #note despite file name these outputs do contain chrZ
  bcftools view -i 'INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o ${OUTNAME}/${chr}_${OUTNAME}_autosomes_dap.vcf.gz $VCF_DIR/${chr}_phased_merge.vcf.gz
  bcftools annotate --set-id '%CHROM\:%POS\:%REF%FIRST_ALT' -O z -o ${OUTNAME}/${chr}_${OUTNAME}_autosomes_dap_ID.vcf.gz ${OUTNAME}/${chr}_${OUTNAME}_autosomes_dap.vcf.gz
  qctool -g ${OUTNAME}/${chr}_${OUTNAME}_autosomes_dap_ID.vcf.gz -ofiletype binary_ped -og ${OUTNAME}/${chr}_${OUTNAME}
  #must enforce numeric chr names
  mv ${OUTNAME}/${chr}_${OUTNAME}.bim ${OUTNAME}/${chr}_${OUTNAME}.orig

  #this doesnt work, has to be actual numbers <95
  ##paste <(cut -f 1 -d " " ${OUTNAME}/${chr}_${OUTNAME}.orig | sed 's/A/99/g' | sed 's/Z/99/g' | sed 's/LGE22/chr29/g' | sed 's/[^0-9]//g') <(cut -d " " -f 2,3,4,5,6 ${OUTNAME}/${chr}_${OUTNAME}.orig) > ${OUTNAME}/${chr}_${OUTNAME}.bim
  #make fake numbers each chr that has letters
  paste <(cut -f 1 -d " " ${OUTNAME}/${chr}_${OUTNAME}.orig | sed 's/1A/31/g' | sed 's/4A/32/g' | sed 's/Z/33/g' | sed 's/LGE22/chr29/g' | sed 's/[^0-9]//g') <(cut -d " " -f 2,3,4,5,6 ${OUTNAME}/${chr}_${OUTNAME}.orig) > ${OUTNAME}/${chr}_${OUTNAME}.bim

  echo ${OUTNAME}/${chr}_${OUTNAME} >> ${OUTNAME}/path_to_all_bfiles.txt
done < all_chr_list.txt
