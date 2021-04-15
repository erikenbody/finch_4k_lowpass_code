
#for every sample, run one job per chromosome for whatshap
#script found in the tool directory

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/23_shapeit4_phasing_new_samps
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/07_BIALLELIC_ONLY/AC_FILTERED/finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC_ACf.vcf.gz

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p $SAMPLE
  cd $SAMPLE

  RES=$(sbatch --array=1-31 -M snowy ~/fc/finch_code/lowpass/23_shapeit4_phasing_new_samps/tools/submit_whatshap_streamlined.sh $SAMPLE)
  echo $RES >> whatshap_job_ids.txt
  sleep 2

  cd $WORK_D

done < vcf_sample_list.txt

#stopped at 100Gen7679 (too many job submissions), will have to come back

#later, decided too many arrays just run it as a loop:

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p $SAMPLE
  cd $SAMPLE

  RES=$(sbatch -M snowy ~/fc/finch_code/lowpass/23_shapeit4_phasing_new_samps/tools/submit_whatshap_streamlined.sh $SAMPLE)
  echo $RES >> whatshap_job_ids.txt
  sleep 2

  cd $WORK_D

done < remaining_samps.txt

#after, see if there right number of VCFs are there

while read SAMPLE
do
  COUNT=$(ls -1 $SAMPLE/*vcf.gz | wc -l)
  echo $SAMPLE $COUNT >> vcf_counts.txt
done < vcf_sample_list.txt
