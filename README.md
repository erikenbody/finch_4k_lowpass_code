# finch_analysis_curated
A repository of some scripts used to process and analyze low coverage whole-genome sequencing data. The main goal of this is to generate a reference haplotype set using a VCF of individuals sequenced at high coveraged. Then, use this reference panel to impute genotypese on a large panel of individuals sequenced at low (2-3x coverrage).

This is moderately organized, but may be useful for others interested in this approach to lowpass genotyping.

The first folder sets up a haplotype reference panel using high coverage sequenced individuals. This requirers a joint-genotyped VCF file, the corresponding BAM files, and the reference genome.
- read backed phasing: https://github.com/whatshap/whatshap
- genotype imputation and phasing: https://github.com/odelaneau/shapeit4
- VCF and BCF manipulation: https://github.com/samtools/bcftools

The second folder takes the informative sites in the reference panel and imputes data for a large number of low coveraged sequenced samples.
- genotype imputation: https://github.com/odelaneau/GLIMPSE

The output is a phased and imputated VCF file for all low coverage individuals
