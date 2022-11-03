library(GenomicRanges)

chr_order <- c("chr1", "chr1A" , "chr2", "chr3", "chr4","chr4A", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11","chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
               "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chrZ",  "chrunknown")

cam.df.an <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes.csv")

input.var <- "p_wald"

gemma.order <- function(gemma.path, input.var){
  df.gemma <- read_tsv(gemma.path, col_names = T, show_col_types = FALSE)
  
  #if chrZ than read it in
  zpat <- gsub("gemma_out/autosomes", "gemma_out/chrZ", gemma.path)
  
  if(file.exists(zpat) == TRUE){
    chrz <- read_tsv(zpat, col_names = T, show_col_types = FALSE)
    chrz$chr <- gsub(-9, "chrZ", chrz$chr)
    chrz <- chrz %>% separate(rs, into = c(NA, NA, "ps"), sep = ":") %>% 
      mutate(ps = as.numeric(ps))
    
    df.gemma <- bind_rows(df.gemma, chrz)
  }
  
  df.gemma$chr_ordered <- factor(df.gemma$chr, levels = chr_order)
  df.gemma <- df.gemma %>% dplyr::arrange(chr_ordered, ps)
  df.gemma$row<-1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  df.gemma$log_p <- -log10(df.gemma[,input.var])
  df.gemma$df.log_p_rollmean <- zoo::rollmean(df.gemma[,input.var],50,fill=NA)
  df.gemma$chr_labels <- gsub("chr", "", df.gemma$chr_ordered)
  chr_breaks <- df.gemma %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row), .groups = "keep")
  df.gemma
}

gcta.order <- function(gemma.path, input.var){
  df.gemma <- read_tsv(gemma.path, col_names = T)
  names(df.gemma) <- gsub("Chr", "chr",
                          gsub("CHR", "chr",
                               gsub("bp", "ps",
                                    gsub("POS", "ps",names(df.gemma)))))
  
  df.gemma$chr <- gsub("31", "1A",
                       gsub("32", "4A",
                            gsub("33", "Z", df.gemma$chr)))
  df.gemma$chr <- paste0("chr",df.gemma$chr)
  
  df.gemma$chr_ordered <- factor(df.gemma$chr, levels = chr_order)
  df.gemma <- df.gemma %>% dplyr::arrange(chr_ordered, ps)
  df.gemma$row<-1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  df.gemma$log_p <- -log10(df.gemma[,input.var])
  df.gemma$df.log_p_rollmean <- zoo::rollmean(df.gemma[,input.var],50,fill=NA)
  df.gemma$chr_labels <- gsub("chr", "", df.gemma$chr_ordered)
  chr_breaks <- df.gemma %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row), .groups = "keep")
  df.gemma
}


peak_list_permutation <- function(input.df, input.var, thresholdA, thresholdZ){
  input.df <- as.data.frame(input.df)
  
  auto_fst <- input.df %>% filter(chr_ordered != "chrZ" & chr_ordered!="chrunknown")
  Z_fst <- input.df %>% filter(chr_ordered == "chrZ")
  
  peaks_fst_auto <- auto_fst %>% filter(!!sym(input.var) > thresholdA)
  
  
  peaks_fst_Z <- Z_fst %>% filter(!!sym(input.var) > thresholdZ)
  peaks_fst <- rbind(peaks_fst_auto, peaks_fst_Z)
  
  #fst_range <- GRanges(peaks_fst$chr, IRanges(as.numeric(peaks_fst$ps), as.numeric(peaks_fst$ps)))
  
  #create granges object for individual SNPs that are above the threshold
  range.peaks_fst <- GRanges(peaks_fst$chr, IRanges(as.numeric(peaks_fst$ps), as.numeric(peaks_fst$ps)))
  
  #merge windows by some distance to other SNPs. defining it here by 75kb
  reduce.range.peaks_fst <- GenomicRanges::reduce(range.peaks_fst, min.gapwidth = 75000) #merge closest window
  #remove ranges that are made up of only 1 SNP. 
  reduce.range.peaks_fst <- reduce.range.peaks_fst[width(reduce.range.peaks_fst) > 1]
  #give each "peak" a name
  names(reduce.range.peaks_fst) <- 1:length(reduce.range.peaks_fst)
  
  datalist <- list()
  #create a dataframe where each SNP gets an peak ID
  for (i in (1:length(reduce.range.peaks_fst))){
    x<-peaks_fst[queryHits(findOverlaps(range.peaks_fst, reduce.range.peaks_fst[i, ])), ]
    x$peak <- i
    datalist[[i]] <- x # add it to your list
  }
  peaks.named = do.call(rbind, datalist)
  
  #an alternative way to make the bed file directly from granges, but it comes with 0 metadata
  #peak.bed <- as.data.frame(reduce.range.peaks_fst)
  #peak.bed$peak <- names(reduce.range.peaks_fst)
  
  return(peaks.named)
  
}

gemma.peak.bed <- function(peaks.named){
  #create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% group_by(chr,peak) %>% 
    summarise(start = min(ps),
              end = max(ps),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
              .groups = "keep")
  return(peak.bed)
}

gemma.peak.ROW.bed <- function(peaks.named){
  #create a bed file that has a summary of each of these regions
  peak.bed <- peaks.named %>% group_by(chr,peak) %>% 
    summarise(start = min(row),
              end = max(row),
              length = end - start,
              num.snps = n(),
              mean.log_p = mean(log_p),
              max.log_p = max(log_p),
              .groups = "keep")
  return(peak.bed)
}

parv.df.an <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes.csv") 
parv.df.exons <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_exons.csv") 

gemma.get_genes <- function(comp.bed){
  parv.df.an.gr <- GRanges(parv.df.an$seqid, IRanges(as.numeric(parv.df.an$start), as.numeric(parv.df.an$end)))
  
  comp.bed <- as.data.frame(comp.bed)
  comp.bed.gr <- GRanges(comp.bed$chr, IRanges(comp.bed$start, comp.bed$end), peak = comp.bed$peak)
  
  comp_overlap <- findOverlaps(comp.bed.gr, parv.df.an.gr)
  
  comp.annotated <- cbind(comp.bed[queryHits(comp_overlap), ], parv.df.an[subjectHits(comp_overlap), c(-1, -2, -3)])
  
  comp.annotated <- comp.annotated %>% mutate(gene_name = gene_symbol) %>% 
    select(-gene_symbol, -Name)
  
  comp.annotated
}

gemma.get_genes.wide <- function(comp.bed){
  parv.df.an.gr <- GRanges(parv.df.an$seqid, IRanges(as.numeric(parv.df.an$start - 1000000), as.numeric(parv.df.an$end + 1000000)))
  
  comp.bed <- as.data.frame(comp.bed)
  comp.bed.gr <- GRanges(comp.bed$chr, IRanges(comp.bed$start, comp.bed$end), peak = comp.bed$peak)
  
  comp_overlap <- findOverlaps(comp.bed.gr, parv.df.an.gr)
  
  comp.annotated <- cbind(comp.bed[queryHits(comp_overlap), ], parv.df.an[subjectHits(comp_overlap), c(-1, -2, -3)])
  
  comp.annotated <- comp.annotated %>% mutate(gene_name = gene_symbol) %>% 
    select(-gene_symbol, -Name)
  
  comp.annotated
}


manc2 <- function(df.in, input.var){
  chr_breaks <- df.in %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row), .groups = "keep")
  
  #for labels, remove > chr20 as these dont plot well
  chr_breaks <- chr_breaks %>% filter(!chr_labels %in% c(21,22,23,24,25,26,27,28))
  
  chrom.colors <- data.frame(chr_ordered = grep("chr", unique(df.in$chr_ordered), value = T),
                             color.num = rep(1:2,length(grep("chr", unique(df.in$chr_ordered))))) %>% 
    distinct(chr_ordered, .keep_all = T)
  
  df.in2 <- df.in %>% #mutate(row = 1:n()) %>% 
    left_join(chrom.colors, by = "chr_ordered") %>% 
    mutate(color.num = as.factor(color.num))
  
  df.in2 %>% 
    filter(chr_labels != "unknown" & !is.na(row)) %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num")) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
          axis.title.x=element_blank(),
          #axis.text.x = element_text(angle = 45, color = "black"),
          #axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=10),
          axis.text = element_text(size=10),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c("grey30","grey70"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    labs(y=input.var) 
  #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
  #                   labels = chr_breaks$chr_labels)
}

#manc3 is for colors
manc3 <- function(df.in, input.var, color1 = "grey30", color2 = "grey70"){
  chr_breaks <- df.in %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row), .groups = "keep")
  
  #for labels, remove > chr20 as these dont plot well
  chr_breaks <- chr_breaks %>% filter(!chr_labels %in% c(21,22,23,24,25,26,27,28))
  
  chrom.colors <- data.frame(chr_ordered = grep("chr", unique(df.in$chr_ordered), value = T),
                             color.num = rep(1:2,length(grep("chr", unique(df.in$chr_ordered))))) %>% 
    distinct(chr_ordered, .keep_all = T)
  
  df.in2 <- df.in %>% #mutate(row = 1:n()) %>% 
    left_join(chrom.colors, by = "chr_ordered") %>% 
    mutate(color.num = as.factor(color.num))

  df.in2 %>% 
    filter(chr_labels != "unknown" & !is.na(row)) %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num")) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
          axis.title.x=element_blank(),
          #axis.text.x = element_text(angle = 45, color = "black"),
          #axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=10),
          axis.text = element_text(size=10),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c(color1,color2))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    labs(y=input.var) 
  #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
  #                   labels = chr_breaks$chr_labels)
}

