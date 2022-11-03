#!/bin/bash
#SBATCH -A snic2022-5-241
#SBATCH -p core -n 1
#SBATCH -M rackham
#SBATCH -t 01-00:00:00
#SBATCH -J gemma
#SBATCH -e gemma_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gemma_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.14 BEDTools/2.29.2

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/PEAK_REGIONS
cd $WORK_D

PEAKS=autosomes_Daphne_cluster1_multivariate_PCsp_lmm4_PEAKS_1mb_ext.bed

GZVCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/tsinfer/autosomes_all_highcov_lowcov_phased_AA.vcf.gz

while read chr start end locus
do
    awk -v num=$locus 'NR == num' $PEAKS > tmp_peak.txt #franken code to ensure that tmp file has tab seperation
    bcftools view -R tmp_peak.txt -O z -o peak_${locus}.vcf.gz $GZVCF
    rm tmp_peak.txt
done < $PEAKS

OUTNAME=Daphne_cluster1
#OUTNAME=Daphne_cluster2
#OUTNAME=Daphne_cluster3
echo $OUTNAME

SUBSET=~/fc/finch_code/lowpass/15_GEMMA/file_input/${OUTNAME}_samples.txt
#for some reason col 2 is species
cut -f 1 $SUBSET > ${OUTNAME}_samples_only.txt

while read chr start end locus
do
    awk -v num=$locus 'NR == num' $PEAKS > tmp_peak.txt #franken code to ensure that tmp file has tab seperation
    bcftools view -S ${OUTNAME}_samples_only.txt -R tmp_peak.txt -O z -o peak_${locus}_${OUTNAME}.vcf.gz $GZVCF
    rm tmp_peak.txt
done < $PEAKS
