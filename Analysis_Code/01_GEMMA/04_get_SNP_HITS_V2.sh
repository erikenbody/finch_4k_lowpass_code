
WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/GEMMA/SNP_HITS
cd $WORK_D

ml bcftools/1.12

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/merged_vcfs/lowcoverage_highcoverage_merge/autosomes_all_highcov_lowcov_phased.vcf.gz

for TOP_SNPS in autosomes_Daphne_cluster1_multivariate_INTs_lmm1_OUTLIER_SNPs_4bcftools.txt #autosomes_Daphne_cluster1_multivariate_INTs_lmm1_OUTLIER_SNPs_TOP100.txt #autosomes_for_ful_mag_species_lmm_OUTLIER_SNPs_TOP100.txt #autosomes_Daphne_cluster1_multivariate_lmm4_OUTLIER_SNPs_TOP100.txt
do
  #TEST_DIR=${TOP_SNPS/_OUTLIER_SNPs_TOP100.txt/}
  TEST_DIR=${TOP_SNPS/_OUTLIER_SNPs_4bcftools.txt/}
  mkdir -p $TEST_DIR
  mkdir -p $TEST_DIR/peak_genotypes
  cut -f 1,2 $TOP_SNPS > $TEST_DIR/peak_genotypes/${TEST_DIR}_snps.txt
  bcftools view -R $TEST_DIR/peak_genotypes/${TEST_DIR}_snps.txt $VCF -O z -o $TEST_DIR/peak_genotypes/${TEST_DIR}_lmm_OUTLIER_SNPs_TOP100.vcf.gz


  VCF=$TEST_DIR/peak_genotypes/${TEST_DIR}_lmm_OUTLIER_SNPs_TOP100.vcf.gz
  SAMPLE=$TEST_DIR
  bcftools query -l $VCF | tr "\n" "\t" > $TEST_DIR/peak_genotypes/${SAMPLE}_header.txt
  echo "" >> $TEST_DIR/peak_genotypes/${SAMPLE}_header.txt #add new line at end of header
  bcftools query -f '[%GT\t]\n' $VCF > $TEST_DIR/peak_genotypes/${SAMPLE}_raw.genos
  cat $TEST_DIR/peak_genotypes/${SAMPLE}_header.txt $TEST_DIR/peak_genotypes/${SAMPLE}_raw.genos > $TEST_DIR/peak_genotypes/${SAMPLE}.genos
  bcftools query -f '%CHROM\t %POS\t %REF\t %ALT\n' $VCF > $TEST_DIR/peak_genotypes/${SAMPLE}_raw.pos
  echo -e "CHROM\tPOS\tREF\tALT" | cat - $TEST_DIR/peak_genotypes/${SAMPLE}_raw.pos > $TEST_DIR/peak_genotypes/${SAMPLE}.pos
  paste -d "\t" $TEST_DIR/peak_genotypes/${SAMPLE}.pos $TEST_DIR/peak_genotypes/${SAMPLE}.genos > $TEST_DIR/peak_genotypes/${SAMPLE}.pos.genos

  sed $TEST_DIR/peak_genotypes/${SAMPLE}.pos.genos \
    -e "s:0|1:1:g" \
    -e "s:1|0:1:g" \
    -e "s:0|0:0:g" \
    -e "s:1|1:2:g" > $TEST_DIR/peak_genotypes/${SAMPLE}.bin.pos.genos

  rm $TEST_DIR/peak_genotypes/${SAMPLE}.genos
  rm $TEST_DIR/peak_genotypes/${SAMPLE}.pos
  rm $TEST_DIR/peak_genotypes/${SAMPLE}_header.txt
  rm $TEST_DIR/peak_genotypes/*_raw*
done
