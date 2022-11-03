library(data.table)
library(tidyverse)
library(patchwork)

#####NOTE CHANGED TO LOG10
#####
fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

source("code/finch_code/lowpass/15_GEMMA/gemma_functions_2021.R")
locus <- "gwas_genotype_chr1A_9"
out.name<-"autosomes_Daphne_cluster1_multivariate_INTs_lmm1"
span <- 500000
  
#get LD structure
lds_seg = read.table("data/gcta/greml_ldms.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
lds_seg$chr <- gsub("31", "1A",
                    gsub("32", "4A",
                         gsub("33", "Z", lds_seg$chr)))
lds_seg$chr <- paste0("chr",lds_seg$chr)
lds_seg$index <- paste(lds_seg$chr, lds_seg$bp, sep = "_")

df_vep <- read.table("data/ensembl_variant_effects/Daphne_cluster1_multivariate_PCsp_significant_snps_VEP.txt")

df_vep <- df_vep %>% 
  #filter(V4!="intergenic_variant") %>% 
  filter(V4=="missense_variant") %>% 
  separate(V2, into = c(NA, "ps"), sep = "-") %>% 
  separate(V1, into = c("chr", NA, NA), sep = "_", remove = F) %>% 
  mutate(ps = as.numeric(ps)) %>% 
  mutate(Consequence = as.factor(V4))

#Create a custom color scale
library(RColorBrewer)
myColors <- rep(brewer.pal(9,"Set1"),2)
names(myColors) <- levels(df_vep$Consequence)
colScale <- scale_colour_manual(name = "grp",values = myColors)

#V2
gwas_loci.og <- c("gwas_genotype_chr1_2_ext", "gwas_genotype_chr1A_17", "gwas_genotype_chr1A_14",
               "gwas_genotype_chr1A_9", "gwas_genotype_chr2_18_ext", 
               "gwas_genotype_chr9_20")


ld_filt_loc.names <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")


pubplot <- function(locus, out.name, span = 500000, trimA = NA, trimB = NA){
  chr.idx <- strsplit(gsub("gwas_genotype_", "", locus), split = "_")[[1]][1]
  peak.id <- strsplit(gsub("gwas_genotype_", "", locus), split = "_")[[1]][2]
  tidy_name <- ld_filt_loc.names[match(locus, ld_filt_loc.names$bed.locus), "tidy_locus"]
  
  df.gemma <- read.table(paste0("output/GEMMA/processed/zoomed_windows/",out.name,"_", peak.id, "_", chr.idx, "_processed.txt"), header = T)
  
  df.gemma.bed.f <- read.table(paste0("output/GEMMA/processed/",out.name,"_PEAKS.bed"))
  names(df.gemma.bed.f) <- c("chr","start","end", "peak")
  df.gemma.bed.f.genes.wide <- gemma.get_genes.wide(df.gemma.bed.f)
  
  #output genes in peaks
  df.gemma.bed.f.genes <-read.table(paste0("output/GEMMA/processed/",out.name,"_lmm_GENES.txt"), header = T)
  
  maxlogp <- max(df.gemma$log_p, na.rm = T)

  print(peak.id)

  #regions
  if(is.na(trimA)){
    df.peak.focus <- df.gemma.bed.f %>% 
      filter(peak == peak.id) %>% 
      mutate(start.ext = start - span, end.ext = end + span)
  } else{
    df.peak.focus <- df.gemma.bed.f %>% 
      filter(peak == peak.id) %>% 
      mutate(start.ext = start + trimA, end.ext = end - trimB,
             start = start.ext + trimA/8, end = end.ext - trimB/8)
  }

  #pvals
  df.gemma.new <- df.gemma %>% 
    filter(chr == df.peak.focus$chr & ps > df.peak.focus$start.ext & ps < df.peak.focus$end.ext)
  #genes
  df.genes.idx <- df.gemma.bed.f.genes.wide #%>% 

    
  df.genes.focus <- left_join(df.genes.idx[,c("gene_id", "peak")], parv.df.an, by = "gene_id") %>% 
    dplyr::rename(CHROM = seqid) %>% 
    mutate(mid = (start + end)/2) %>% 
    filter(peak == peak.id)
  df.genes.focus$gene_symbol <- toupper(df.genes.focus$gene_symbol)
  
  df.genes.focus <- df.genes.focus %>% 
    filter(CHROM == df.peak.focus$chr & start > (df.peak.focus$start - 50000) & end < (df.peak.focus$end + 50000)) 

  #for arrow direction
  df.genes.focus <- df.genes.focus %>% 
    mutate(start.strand = ifelse(strand == "-", end, start),
           end.strand = ifelse(strand == "-", start - 20000, end + 20000))
  
  #recombination
  #df_recomb.focus <- df_recomb %>% 
  #  filter(chr == chr.idx) %>% 
  #  filter(ps > df.peak.focus$start.ext & ps < df.peak.focus$end.ext)
  
  #ld
  lds_seg.focus <- lds_seg %>% filter(chr == chr.idx & bp > df.peak.focus$start.ext & bp < df.peak.focus$end.ext)
  
  df.gemma.new <- left_join(df.gemma.new, df_vep[,c("chr", "ps", "Consequence")], by = c("chr", "ps"))
  
  #get exons
  df.exons.focus <- left_join(df.genes.focus[,c("gene_id", "peak")], parv.df.exons, by = "gene_id")

  ##plotting
  
  #smart tip for aligning labels here: https://stackoverflow.com/questions/51024675/preserving-order-with-geom-text-repel
  #I try to add a bit of a buffer so that more can be plotted
  
  
  df.genes.focus$i <- seq(min(df.genes.focus$mid) - 500000, max(df.genes.focus$mid) + 500000, length.out = nrow(df.genes.focus))
  
  df.genes.focus <- df.genes.focus %>% 
    mutate(gene_symbol = ifelse(grepl("ENSC", gene_symbol), "UNK", gene_symbol))
  
  null_theme <- theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), 
                      #axis.line = element_line(colour = "black"),
                      axis.line = element_blank(),
                      axis.title.x=element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.title.y=element_blank(), 
                      axis.text.y = element_blank(), 
                      axis.ticks.y = element_blank(),
                      legend.position = "none")
  
  p.genes <- ggplot(data = df.gemma.new, aes(x=ps/10000,y=value)) +
    geom_segment(data = df.exons.focus, aes(y = 2, yend = 2, x = start/1000000, xend = end/1000000, color = strand), 
                 size = 5) +
    geom_segment(data = df.genes.focus, aes(y = 2, yend = 2, x = (start.strand - 50)/1000000, xend = (end.strand + 50)/1000000, color = strand), 
                 size = 1, linejoin='mitre', arrow = arrow(length = unit(0.1, "cm"))) +
    geom_text(
      data = df.genes.focus,
      mapping = aes(y = 2.1, x = i/1000000, label = gene_symbol),
      parse = TRUE, hjust = 0, angle = 90, size = 5
    ) +
    geom_segment(
      data = df.genes.focus,
      mapping = aes(y = 2, yend = 2.1, x = mid/1000000, xend = i/1000000),
      size = 0.1
    ) +
    null_theme +
    ylim(1.9,2.7) +
    xlim(min(df.gemma.new$ps, na.rm = T)/1000000, max(df.gemma.new$ps,na.rm = T)/1000000) +
    scale_color_manual(values = c("red","blue"))
  
  p.genes
  
  p.gemma.outlier <- df.gemma.new %>%
    filter(log_p > 7.7)
  
  
  
  df.conseq <- df.gemma.new %>% filter(!is.na(Consequence))
  df.conseq$Consequence <- droplevels(as.factor(df.conseq$Consequence))
  
  #commented out the missense mutations for now
  p.gemma <- ggplot(data = df.gemma.new, aes(x = ps/1000000, y = log_p)) + 
    geom_point(data = df.gemma.new, size = 1) +
    #geom_point(data = df.conseq, aes(color = Consequence), size = 1.5) +
    geom_hline(yintercept = 7.7, linetype =2, color = "red") + 
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    theme(axis.title.x=element_blank(),
          text = element_text(size = 18))+
    labs(y=expression(paste("-log"[10], italic("P"),"-value")), x = unique(df.peak.focus$chr)) +
    ggtitle(paste0(unique(df.peak.focus$chr),", ", "locus = ",tidy_name)) +
    ylim(0,maxlogp) +
    colScale
  
  p.gemma <- p.gemma + theme(legend.position = "none")
  p.gemma
  
  lds_seq.chr <- lds_seg %>%
    filter(chr == chr.idx) %>% 
    mutate(index = 1:n()) 
  
  #remake peak list from original bed
  df.peak.focus2 <- df.gemma.bed.f %>% 
    filter(peak == peak.id) 
  
  axis1 <- lds_seq.chr[match(df.peak.focus2$start, unique(lds_seq.chr$bp)), "index"]
  axis2 <- lds_seq.chr[match(df.peak.focus2$end, unique(lds_seq.chr$bp)), "index"]
  
  #this really didnt want to work
  #have to have x as a true numeric (not the 'ps' with missign in between values)
  #also raster not tile for saving (otherwise colors are wrong in output)
  
  p.ldchr <- lds_seq.chr %>%
    ggplot(aes(x = index, y = chr, fill = ldscore_region)) + 
    #geom_tile() +
    geom_raster()+
    geom_vline(xintercept = axis1) +
    geom_vline(xintercept = axis2) +
    scale_fill_gradient(low = "lightblue", high = "red") +
    null_theme 
    
  
  layout <- "
  A
  A
  A
  B
  B
  B
  B
  C
  "
  
  #pdf(paste0("output/GEMMA/processed/zoomed_windows/",out.name,peak.id, "_", chr.idx, "_zoomed.pdf"), width = 8.5, height = 10)
  pcomb <- p.genes / p.gemma / p.ldchr  +
    plot_layout(design = layout)
  #dev.off()
  pcomb
  ggsave(paste0("output/GEMMA/processed/zoomed_windows/",out.name,peak.id, "_", chr.idx, "_zoomed.pdf"), width = 8.5, height = 6, dpi = 320)
  ggsave(paste0("output/GEMMA/processed/zoomed_windows/",out.name,peak.id, "_", chr.idx, "_zoomed.png"), width = 10, height = 6, dpi = 320)
  
  return(pcomb)
}
 
for (locus in gwas_loci.og){
  pubplot(locus, "autosomes_Daphne_cluster1_multivariate_PCsp_lmm4", 500000, NA, NA)
}

pubplot("gwas_genotype_chr1A_9", "autosomes_Daphne_cluster1_multivariate_PCsp_lmm4", 500000, 1.8e6, 0.9e6)


