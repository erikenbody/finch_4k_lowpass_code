
#now this working D is where run with recombination map will happen
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse
cd $WORK_D

#11 dec TESTING
#sbatch --array=1-25 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh chr3
#checked these results and memory capped at 18gb, but it successfully used all 10 threads

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-25 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh $CHROM

  sleep 2

  cd $WORK_D

done < chr_list.txt

#removed chr3 (already ran) and moved chr1A to start of chr_list
#smaller chromosomes should run faster, so running these with 5 cores instead:

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-25 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools_small_chrom.sh $CHROM

  sleep 2

  cd $WORK_D

done < smaller_chr_list.txt


#####SUBMIT 2020 SAMPLES#####
##All I did here is replace the location of the interval files (from bcftools) in the run_glimpse script
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/lowpass_round2_data
cd $WORK_D

##sbatch --dependency=afterok:17293338 --array=1-12 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh chr3

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-12 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh $CHROM

  sleep 2

  cd $WORK_D

done < /crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/large_chr_list.txt

#smaller chromosomes should run faster, so running these with 5 cores instead:

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-12 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools_small_chrom.sh $CHROM

  sleep 2

  cd $WORK_D

done < /crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/smaller_chr_list.txt

########February 2021 samples, last round########

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/lowpass_round3_data
cd $WORK_D

##test
##sbatch -p devel -t 00:10:00 --array=1 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh chr3

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-6 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools.sh $CHROM

  sleep 2

  cd $WORK_D

done < /crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/large_chr_list.txt

#smaller chromosomes should run faster, so running these with 5 cores instead:

while read CHROM
do
  echo "$CHROM"
  mkdir -p $CHROM
  cd $CHROM

  sbatch --array=1-6 ~/fc/finch_code/lowpass/14_external_imputation_glimpse/tools/run_glimpse_bcftools_small_chrom.sh $CHROM

  sleep 2

  cd $WORK_D

done < /crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/smaller_chr_list.txt

##I reccomend in the future to set this up in some way that records the job id when submitting, so they can easily be canceled
