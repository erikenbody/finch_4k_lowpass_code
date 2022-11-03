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

plot.inflation <- function (x, size = 2) {

  # Get the number of p-values.
  n <- length(x)
  
  # Compute the negative log10(p-values), and sort them from largest
  # to smallest.
  y <- rev(sort(-log10(x)))
  df.pvals <- data.frame(x = -log10((1:n)/n),y = y)
  
  z <- tail(y, n = 0.9*(length(y)))
  n2 <- length(z)
  df.f90 <- data.frame(x = -log10((1:n2)/n2),y = z)
  lm.df90 <- lm(y ~ x, data = df.f90)
  df.f90$fitted.y <- lm.df90$fitted.values
  
  qq.plot <- ggplot() +
    geom_line(data = df.f90, aes(x = x, y = fitted.y)) +
    geom_abline(intercept = 0,slope = 1,color = "magenta") +
    geom_point(data = df.pvals, aes(x = x, y = y), color = "dodgerblue",shape = 20,size = 2) +
    labs(x = "Expected -log10 p-value",
         y = "Observed -log10 p-value") +
    theme(axis.line = element_blank())
  
  return(qq.plot)

}

# loop over runs ---------------# loop over runs ---------------# loop over runs ---------------
all_cluster.l <- grep("backup",list.files(path = "data/GEMMA/", pattern = "assoc.txt", full.names = T, recursive = T), value = T, invert = T)
names_cluster.l <- grep("chr", grep("lmm", all_cluster.l, value = T), invert = T, value = T)

run_name <- "data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_multivariate_PCsp_lmm4.assoc.txt"

process_gemma <- function(run_name, threshold, threshold.strict, gwas_stat, color1, color2){
  print(run_name)
  run_name2 <- basename(run_name)
  out.name <- gsub(".assoc.txt","", run_name2)
  
  gemma_files <- grep(run_name, all_cluster.l, value = T)
  
  df.gemma <- gemma.order(gemma_files, gwas_stat)
  write.table(df.gemma, paste0("output/GEMMA/processed/",out.name,"_processed.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  #remove low log p to reduce memory footprint. 
  df.gemma.filt <- df.gemma %>% filter(log_p > 1)
  ##Also remove MAF < 0.05. negligible difference to plotting
  #df.gemma.filt <- df.gemma %>% filter(log_p > 1 & af > 0.05 & af < 0.95)
  
  #define peaks using threshold set by permutation test (external gemma_wrapper script)
  #internally, this merges windows within 75kb of each other and removes single snp regions
  df.gemma.peaks <- peak_list_permutation(df.gemma.filt, "log_p",threshold,threshold)
  #define intervals where the peaks are, keeping information on length and mean/max p values
  df.gemma.bed <- gemma.peak.bed(df.gemma.peaks)
  #use a genome-wide index ("row") for manhattan plotting
  df.gemma.ROW.bed <- gemma.peak.ROW.bed(df.gemma.peaks)
  
  #create a filtered bed file by some custom parameters 
  df.gemma.bed.f <- df.gemma.bed %>% 
    filter(max.log_p > threshold.strict & length > 100000) 
  
  df.gemma.ROW.bed.f <- df.gemma.ROW.bed %>% 
    filter(peak %in% df.gemma.bed.f$peak)
  
  if(nrow(df.gemma.bed.f) > 0){
  df.gemma.bed.f <- df.gemma.bed.f%>% 
    ungroup() %>% 
    mutate(chr_ordered = factor(chr, levels = chr_order)) %>% #arrange by chr order
    arrange(chr_ordered, start) %>% #arrange by chr order
    mutate(peak.original = peak, #new col for original peak name
           peak.new = 1:n(), #new names for peaks
           peak = peak.new) #reset original peak coloumn to the renumbered version
  
  df.gemma.ROW.bed.f <- df.gemma.ROW.bed.f %>% 
    ungroup() %>% 
    mutate(chr_ordered = factor(chr, levels = chr_order)) %>% #arrange by chr order
    arrange(chr_ordered, start) %>% #arrange by chr order
    mutate(peak.original = peak, #new col for original peak name
           peak.new = 1:n(), #new names for peaks
           peak = peak.new) #reset original peak coloumn to the renumbered version
  
  peak.name.simple <- df.gemma.ROW.bed.f %>% select(peak.new, peak.original)
  
  }
  
  
  #get genes
  df.gemma.bed.f.genes <- gemma.get_genes(df.gemma.bed.f)
  df.gemma.bed.f.genes.wide <- gemma.get_genes.wide(df.gemma.bed.f)
  
  #a method to visualize quickly the peaks per chr
  #test.chr <- "chr1A"
  #input.df %>% filter(chr == test.chr & log_p > 4) %>% 
  #  ggplot() + geom_point(aes(x = ps, y = log_p)) +
  #  geom_segment(data = base::subset(df.gemma.bed.f, chr %in% test.chr), aes(x = start, xend = end, y = 30, yend = 30), color = "blue", size = 2)
  
  #output datafiles
  #create a SNP file that is only the peaks after filtering
  df.gemma.peaks.out <- df.gemma.peaks %>% filter(peak %in% df.gemma.bed.f$peak.original)
  df.gemma.peaks.out <- left_join(df.gemma.peaks.out, peak.name.simple, by = c("peak" = "peak.original")) %>% 
    dplyr::rename(peak.original = peak, peak= peak.new)
  write.table(df.gemma.peaks.out,paste0("output/GEMMA/processed/",out.name,"_OUTLIER_SNPs.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  df.gemma.peaks.out2 <- df.gemma.peaks.out %>% 
    group_by(peak) %>% 
    ungroup() %>% 
    mutate(variant = paste0(allele1, "/",allele0)) %>% 
    mutate(ps2 = ps) %>% 
    select(chr, ps, ps2, variant) %>% 
    arrange(chr, ps)
  write.table(df.gemma.peaks.out2,paste0("output/GEMMA/processed/",out.name,"_OUTLIER_SNPs_4ensembl.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  #create a file that is only the top 100 SNPs per peak
  df.gemma.peaks.slice <- df.gemma.peaks.out %>% group_by(peak) %>% 
    slice_max(log_p, n = 100) %>% 
    select(chr, ps, log_p, peak)
  write.table(df.gemma.peaks.slice,paste0("output/GEMMA/processed/",out.name,"_OUTLIER_SNPs_TOP100.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  df.gemma.peaks.slice2 <- df.gemma.peaks.out %>% 
    group_by(peak) %>% 
    slice_max(log_p, n = 100) %>% 
    ungroup() %>% 
    mutate(variant = paste0(allele1, "/",allele0)) %>% 
    mutate(ps2 = ps) %>% 
    select(chr, ps, ps2, variant) %>% 
    arrange(chr, ps)
  write.table(df.gemma.peaks.slice2,paste0("output/GEMMA/processed/",out.name,"_OUTLIER_SNPs_TOP100_4ensembl.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  
  #output bed file of UNFILTERED peaks
  df.gemma.bed %>% select(chr,start,end, peak) %>% 
    write.table(paste0("output/GEMMA/processed/",out.name,"PEAKS_unfiltered.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  #output bed file of peaks
  df.gemma.bed.f %>% select(chr,start,end, peak) %>% 
    write.table(paste0("output/GEMMA/processed/",out.name,"_PEAKS.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  #output genes in peaks
  df.gemma.bed.f.genes %>% 
    write.table(paste0("output/GEMMA/processed/",out.name,"_lmm_GENES.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  maxlogp <- max(df.gemma$log_p, na.rm = T)
  
  
  plot.list <- list()
  pdf(paste0("output/GEMMA/processed/",out.name,"_lmm_per_peak_plot.pdf"), width = 8.5, height = 6, useDingbats = F)
  for (peak.id in unique(df.gemma.bed.f$peak)){
    print(peak.id)
    #regions
    df.peak.focus <- df.gemma.bed.f %>% filter(peak == peak.id) %>% 
      mutate(start.ext = start - 1000000, end.ext = end + 1000000)
    #pvals
    df.gemma.new <- df.gemma %>% filter(chr == df.peak.focus$chr & ps > df.peak.focus$start.ext & ps < df.peak.focus$end.ext)
    #dir.create("output/GEMMA/processed/zoomed_windows")
    write.table(df.gemma.new, paste0("output/GEMMA/processed/zoomed_windows/", out.name, "_", peak.id, "_", unique(df.peak.focus$chr), "_processed.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
    #genes
    df.genes.idx <- df.gemma.bed.f.genes.wide %>% 
      filter(chr == df.peak.focus$chr & start > df.peak.focus$start.ext & end < df.peak.focus$end.ext) %>% 
      select(gene_id, peak)
    df.genes.focus <- left_join(df.genes.idx, parv.df.an, by = "gene_id") %>% 
      dplyr::rename(CHROM = seqid) %>% mutate(mid = (start + end)/2) %>% 
      filter(peak == peak.id)
    df.genes.focus$gene_symbol <- toupper(df.genes.focus$gene_symbol)
    #for arrow direction
    df.genes.focus <- df.genes.focus %>% 
      mutate(start.strand = ifelse(strand == "-", end, start),
             end.strand = ifelse(strand == "-", start - 20000, end + 20000))
    
    #get exons
    df.exons.focus <- left_join(df.genes.idx, parv.df.exons, by = "gene_id")
    
    ##plotting

    
    if (dim(df.genes.focus)[1] == 0) {
      p.genes <- ggplot(data = df.gemma.new, aes(x=ps/10000,y=log_p)) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              #axis.line = element_line(colour = "black"),
              axis.line = element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.y=element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()) +
        ylim(1.9,2.6) +
        xlim(min(df.gemma.new$ps, na.rm = T)/1000000, max(df.gemma.new$ps,na.rm = T)/1000000)
    } else {
      
      #smart tip for aligning labels here: https://stackoverflow.com/questions/51024675/preserving-order-with-geom-text-repel
      #I try to add a bit of a buffer so that more can be plotted
      df.genes.focus$i <- seq(min(df.genes.focus$mid) - 500000, max(df.genes.focus$mid) + 500000, length.out = nrow(df.genes.focus))
      
      p.genes <- ggplot(data = df.gemma.new, aes(x=ps/10000,y=value)) +
        geom_segment(data = df.exons.focus, aes(y = 2, yend = 2, x = start/1000000, xend = end/1000000, color = strand), 
                     size = 5) +
        geom_segment(data = df.genes.focus, aes(y = 2, yend = 2, x = (start.strand - 50)/1000000, xend = (end.strand + 50)/1000000, color = strand), 
                     size = 1, linejoin='mitre', arrow = arrow(length = unit(0.1, "cm"))) +
        geom_text(
          data = df.genes.focus,
          mapping = aes(y = 2.1, x = i/1000000, label = gene_symbol),
          parse = TRUE, hjust = 0, angle = 90, size = 3
        ) +
        geom_segment(
          data = df.genes.focus,
          mapping = aes(y = 2, yend = 2.1, x = mid/1000000, xend = i/1000000),
          size = 0.1
        ) +
        theme(panel.grid.major = element_blank(), 
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
              legend.position = "none") +
        ylim(1.9,2.7) +
        xlim(min(df.gemma.new$ps, na.rm = T)/1000000, max(df.gemma.new$ps,na.rm = T)/1000000) +
        scale_color_manual(values = c("red","blue"))
    }  
    p.genes
    
    p.gemma.outlier <- df.gemma.new %>%
      filter(log_p > 7.7)
    
    p.gemma <- ggplot() + 
      geom_point(data = df.gemma.new, aes(x = ps/1000000, y = log_p), size = .5) +
      geom_point(data = p.gemma.outlier, aes(x = ps/1000000, y = log_p), size = .5, color = "red") +
      geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
      geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
      theme_bw() + 
      #general.theme +
      theme(axis.title.x=element_blank())+
            #axis.text.x=element_blank()) +
      labs(y=expression(paste("-log"[10], italic("P"),"-value")), x = unique(df.peak.focus$chr)) +
      #scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") + 
      ggtitle(paste0(unique(df.peak.focus$chr),", ", "peak = ",peak.id)) +
      ylim(0,maxlogp) 
    
    
    
    
    plot.list[[peak.id]] <- p.genes / p.gemma
    print(plot.list[[peak.id]])
    
  }
  dev.off()
  
  #plot per chr
  #plot.list <- list()
  #pdf(paste0("output/GEMMA/processed/",out.name,"_lmm_per_chr_plot.pdf"), width = 8.5, height = 6, useDingbats = F)
  #for (chr.set in unique(df.gemma.bed.f$chr)){
  #  print(chr.set)
  #  df.bed.new <- df.gemma.bed.f %>% filter(chr == chr.set) %>% 
  #    mutate(start.ext = start - 1000000, end.ext = end + 1000000)
  #  df.gemma.new <- df.gemma.filt %>% filter(chr == chr.set)
  #  
  #  plot.list[[chr.set]] <- ggplot() + 
  #    geom_point(data = df.gemma.new, aes(x = ps, y = log_p)) +
  #    geom_segment(data = df.bed.new, aes(x = start, xend = end, y = 15, yend = 15), color = "blue", size = 2) +
  #    theme_bw() + ggtitle(chr.set) +
  #    geom_hline(yintercept = threshold.strict, color = "red")
  #  print(plot.list[[chr.set]])
  #  
  #}
  #dev.off()
  
  #p.gemma <- manc2(df.gemma.filt, "log_p")
  p.gemma <- manc3(df.gemma.filt, "log_p", color1, color2)
  #visualize peaks whole genome
  #png(paste0("output/GEMMA/processed/",out.name,"_man_peaks.png"), width = 2400, height = 800, res = 200) 
  p.gemma2 <- p.gemma +
    #geom_segment(data = df.gemma.ROW.bed.f, aes(x = start, xend = end, y = max.log_p + 5, yend = max.log_p + 5), color = "blue", size = 2) +
    #geom_hline(yintercept = threshold, linetype = 2, color = "red", size = 0.5) +
    geom_hline(yintercept = threshold.strict, color = "red", size = 0.5) +
    labs(y = expression(paste("-log"[10], italic("P"),"-value")))
  p.gemma2
  ggsave(paste0("output/GEMMA/processed/",out.name,"_man_peaks.png"), width = 10, height = 3, dpi = 300)
  #dev.off()
  
}

#multivariatePC1.sp and PC2.sp
process_gemma("data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_multivariate_PCsp_lmm4.assoc.txt", 6.6, 7.7,"p_wald", fort.color, "#005685")
process_gemma("data/GEMMA//Daphne_cluster2/autosomes_Daphne_cluster2_gemma_out/autosomes_Daphne_cluster2_multivariate_PCsp_lmm4.assoc.txt", 6.6, 7.7,"p_wald", scan.color, "#006346")
process_gemma("data/GEMMA//Daphne_cluster3/autosomes_Daphne_cluster3_gemma_out/autosomes_Daphne_cluster3_multivariate_PCsp_lmm4.assoc.txt", 6.6, 7.7,"p_wald", magn.color, "#906300")

#with hmga2 as covariate for multivariate PCs
process_gemma("data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_PC1_cov.p33100090_lmm4.assoc.txt", 6.6, 7.7,"p_wald", fort.color, "#005685")

#without weight as covariate for multivariate PCs
process_gemma("data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_multivariate_PCsp_no_weight_lmm4.assoc.txt", 6.6, 7.7,"p_wald", fort.color, "#005685")

#body weight with int
process_gemma("data/GEMMA//Daphne_cluster1/autosomes_Daphne_cluster1_gemma_out/autosomes_Daphne_cluster1_weight_INTs_no_beak_cov_lmm4.assoc.txt", 6, 7.7, "p_wald", fort.color, "#005685")
process_gemma("data/GEMMA//Daphne_cluster2/autosomes_Daphne_cluster2_gemma_out/autosomes_Daphne_cluster2_weight_INTs_no_beak_cov_lmm4.assoc.txt", 6, 7.7, "p_wald",  scan.color, "#006346")
process_gemma("data/GEMMA//Daphne_cluster3/autosomes_Daphne_cluster2_gemma_out/autosomes_Daphne_cluster3_weight_INTs_no_beak_cov_lmm4.assoc.txt", 6, 7.7, "p_wald", magn.color, "#906300")
