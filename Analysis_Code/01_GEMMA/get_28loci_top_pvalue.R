library(tidyverse)
library(GenomicRanges)
df.bed28 <- read.table("/Users/erikenbody/Google\ Drive/My\ Drive/Uppsala/Projects/inProgress/DarwinsFinchAssembly/output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS.bed") 
df.fortis.gemm <- read.table("data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_multivariate_INTs_lmm1.assoc.txt", header = T)
df.fortis.gemm$logp <- -log10(df.fortis.gemm$p_wald)
df.gemm.gr <- GRanges(df.fortis.gemm$chr, IRanges(as.numeric(df.fortis.gemm$ps), as.numeric(df.fortis.gemm$ps)), logp = df.fortis.gemm$logp)

comp.bed <- as.data.frame(df.bed28)
comp.bed.gr <- GRanges(comp.bed$V1, IRanges(comp.bed$V2, comp.bed$V3), peak = comp.bed$V4)

comp_overlap <- findOverlaps(comp.bed.gr, df.gemm.gr)

comp.annotated <- cbind(comp.bed[queryHits(comp_overlap), ], df.gemm.gr[subjectHits(comp_overlap)])
head(comp.annotated)

comp.annotated %>% 
  group_by(V4) %>% 
  slice_max(logp, n = 1) %>% 
  #print(n = 28) %>% 
  write_csv("output/GEMMA/asm28loci_top_pvalue.csv")
