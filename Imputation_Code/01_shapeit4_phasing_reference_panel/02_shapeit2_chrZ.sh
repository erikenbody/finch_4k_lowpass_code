#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p node
#SBATCH -C mem256GB
#SBATCH -M rackham #rackham only
#SBATCH -t 03-00:00:00
#SBATCH -J shapeit2
#SBATCH -e shapeit2_%J_%A_%a.err
#SBATCH -o shapeit2_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools bcftools/1.12 plink/1.90b4.9 SHAPEIT/v2.r904

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/23_shapeit4_phasing_new_samps
cd $WORK_D


INTERVAL=chrZ
VCF=${INTERVAL}.415samps.bcf

##Failed attempts to get shapeit 2 to work with --chrX
#bcftools view -O z -o ${INTERVAL}.415samps.vcf.gz $VCF
# bcftools view -i 'F_MISSING < 0.1' -O z -o ${INTERVAL}.415samps_fil2.vcf.gz $VCF
# bcftools view -i 'AC > 5' -O z -o ${INTERVAL}.415samps_fil3.vcf.gz $VCF
# bcftools view -i 'F_MISSING < 0.8' -O z -o ${INTERVAL}.415samps_fil4.vcf.gz $VCF

##Only this was used
# bcftools view -S males_only_samples.txt -O z -o ${INTERVAL}.276samps_males.vcf.gz $VCF
# bcftools view -i 'F_MISSING < 0.8' -O z -o ${INTERVAL}.276samps_males_fil.vcf.gz ${INTERVAL}.276samps_males.vcf.gz

##it doesnt work to manually edit the .fam file to add sex. Must run plink with the sex IDs using update-sex
##shapeit gives error error ERROR: 1 fully missing SNPs (=100%), but only with my custom fam file, not the original fam file...
##Note this file has flipped sex 1 and 2 to reflect that shapeit needs the heterogametic sex to be 1
# plink --vcf ${INTERVAL}.415samps_fil4.vcf.gz --make-bed --out  ${INTERVAL}_fil4.415samps --allow-extra-chr --update-sex finch_with_sex_for_shapeit.txt
# plink --vcf ${INTERVAL}.276samps_males_fil.vcf.gz --make-bed --out  ${INTERVAL}.276samps_males.vcf.gz --allow-extra-chr --update-sex finch_with_sex_for_shapeit.txt


# shapeit -B chrZ.276samps_males.vcf.gz \
#         -M gmap_files/${INTERVAL}_parv.gmap \
#         -O ${INTERVAL}_males_shapeit2.phased \
#         --thread 1 \
#         --force \
#         --chrX

#in the end, I decided to only phase male samples, and leave off females

####Only using this:
module load bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.12
SI=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/shapeit4-4.2.1/bin
#bcftools index ${INTERVAL}.276samps_males.vcf.gz
$SI/shapeit4.2 --input ${INTERVAL}.276samps_males.vcf.gz \
          --output ${INTERVAL}.415samps.postrecomb.bcf \
          --sequencing \
          --map gmap_files/${INTERVAL}_parv.gmap \
          --thread 20 \
          --region $INTERVAL \
          --log ${INTERVAL}.phasing.postrecomb.log
