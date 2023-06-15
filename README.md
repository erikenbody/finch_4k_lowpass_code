[![DOI](https://zenodo.org/badge/358115636.svg)](https://zenodo.org/badge/latestdoi/358115636)

# Imputation code
A repository of some scripts used to process and analyze low coverage whole-genome sequencing data. The main goal of this is to generate a reference haplotype set using a VCF of individuals sequenced at high coveraged. Then, use this reference panel to impute genotypese on a large panel of individuals sequenced at low (2-3x coverrage).

This is moderately organized, but may be useful for others interested in this approach to lowpass genotyping.

The first folder sets up a haplotype reference panel using high coverage sequenced individuals. This requires a joint-genotyped VCF file, the corresponding BAM files, and the reference genome. The software used is:
- read backed phasing: https://github.com/whatshap/whatshap
- genotype imputation and phasing: https://github.com/odelaneau/shapeit4
- VCF and BCF manipulation: https://github.com/samtools/bcftools

The second folder takes the informative sites in the reference panel and imputes data for a large number of low coveraged sequenced samples. The software used is:
- genotype imputation: https://github.com/odelaneau/GLIMPSE

The output is a phased and imputated VCF file for all low coverage individuals

In order to phase samples in the reference panel on chrZ, we only used male samples. These were identified using the snakemake pipeline in sex_ID_using_depth_snakemake. This requires samtools https://github.com/samtools/samtools and compares sequencing depth on chr4 (a moderately sized chromosome) to that on chrZ - producing the output per sample with the ratio of Z to 4 sequencing depth. 

# Other analysis code

Code used as part of manuscript in review "Large effect loci have a prominent role in Darwin’s finch evolution"

bioRxiv preprint:
https://www.biorxiv.org/content/10.1101/2022.10.29.514326v1
