#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 5
#SBATCH -M rackham
#SBATCH -t 2-00:00:00
#SBATCH -J concat_glimpse
#SBATCH -e concat_glimpse_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o concat_glimpse_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs
cd $WORK_D

SUBSET=autosomes

LST=${SUBSET}.list
# ls chr*_phased_merge.vcf.gz | grep -v "chrZ" > ${LST}
#
# bcftools concat -n --threads 5 -f ${LST} -O z -o glimpse_${SUBSET}_phased_merge.vcf.gz
# bcftools index glimpse_${SUBSET}_phased_merge.vcf.gz

# LST=${SUBSET}_ligated.list
# ls chr*_ligated_merge.vcf.gz | grep -v "chrZ" > ${LST}
#
# bcftools concat -n --threads 5 -f ${LST} -O z -o glimpse_${SUBSET}_ligated_merge.vcf.gz
# bcftools index glimpse_${SUBSET}_ligated_merge.vcf.gz

##set up some kind of filtered subset to test our e.g. pca
##VCF=glimpse_${SUBSET}_ligated_merge.vcf.gz
##bcftools query -l $VCF | grep "Dap" > Dap_samples.txt
##bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2 || ALT="*"' -m2 -M2 -S Dap_samples.txt -i 'INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o autosomes_dap_filt.vcf.gz $VCF
##only get biallelics and set some rough AF filters (0.05) and filters that probably are associated with bad phasing, based on examining these distributions
##bcftools view -m2 -M2 -S Dap_samples.txt -i 'AC!=0 & AC!=AN & F_MISSING < 0.2 & ALT!="*" & INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o autosomes_dap_filt.vcf.gz $VCF
##bcftools +prune -w 1000bp -n 1 -O z -o autosomes_dap_filt_1kb_pruned.vcf.gz autosomes_dap_filt.vcf.gz

#just use phase?
#phased
VCF=glimpse_${SUBSET}_phased_merge.vcf.gz
bcftools query -l $VCF | grep "Dap" > Dap_samples.txt
#not using vars calculated on the fly by bcftools
bcftools view -S Dap_samples.txt -i 'INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o autosomes_dap_phased_filt.vcf.gz $VCF
bcftools +prune -w 20000bp -n 1 -O z -o autosomes_dap_phased_filt_pruned_20kb.vcf.gz autosomes_dap_phased_filt.vcf.gz
bcftools index autosomes_dap_phased_filt_pruned_20kb.vcf.gz

#chrZ
bcftools view -S Dap_samples.txt -i 'INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o chrZ_dap_phased_filt.vcf.gz chrZ_phased_merge.vcf.gz
bcftools index chrZ_dap_phased_filt.vcf.gz

#non daphne
bcftools query -l $VCF | grep -v "Dap" > non_Dap_samples.txt
#not using vars calculated on the fly by bcftools
bcftools view -S non_Dap_samples.txt -i 'INFO/RAF>0.05 & INFO/RAF<0.95 & INFO/INFO>0.0625 & INFO/AF>0.05 & INFO/AF<0.95' -O z -o autosomes_NONdap_phased_filt.vcf.gz $VCF
