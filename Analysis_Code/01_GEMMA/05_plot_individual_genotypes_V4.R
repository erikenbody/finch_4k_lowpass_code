#library(data.table)
library(tidyverse)
library(readxl)
#library(dendextend)
library(patchwork)
library(data.table)
library(cluster)
#source("code/finch_code/lowpass/15_GEMMA/gemma_functions_2021.R")

# prefix ------------------------------------------------------------------

prefix <- "autosomes_Daphne_cluster1_multivariate_PCsp_lmm4"

# load data ---------------------------------------------------------------

#load database with lowpass morphometrics
V2.mdb <- read_excel("/Users/erikenbody/Dropbox/Shared_Uppsala/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
V2.mdb <- V2.mdb %>% 
  filter(!is.na(Lowcov_project)) %>% 
  #filter(!is.na(Final.proj) | !is.na(lowpass.2020)) %>% 
  mutate(bill.size = bill.length + bill.width + bill.depth)

##Get PCs and add to V2.md
dap.pcs <- read_tsv("code/finch_code/lowpass/15_GEMMA/file_input/sample_pheno_all_3957_pheno.txt", col_names = F)
dap.pcs.names <- read_tsv("code/finch_code/lowpass/15_GEMMA/file_input/sample_pheno_all_3957_samples.txt", col_names = F)
dap.pcs2 <- dap.pcs[c(19,20)]
dap.pcs2$meaningful.unique <- dap.pcs.names$X1
names(dap.pcs2) <- c("bill.PC1", "bill.PC2", "meaningful.unique")

V2.mdb <- left_join(V2.mdb, dap.pcs2, by = "meaningful.unique")

#load genotypes
df.genos <- read_tsv(paste0("data/GEMMA/SNP_HITS/", prefix, ".bin.pos.genos"))

df.genos <- df.genos[,1:(length(df.genos)-1)]#there is a trailing NA record
#this file looks like this:
#     CHROM      POS REF ALT 01Dap19024 01Dap19030 01Dap19032 01Dap19034 01Dap19035 01Dap19036
#1:  chr1 61973061   T   G          0          0          0          0          0          0
#2: chr1A 30058647   G   A          1          2          0          2          2          2
#3: chr1A 30065178   T   C          1          2          0          2          2          2
#4: chr1A 30324767   T   C          1          0          2          0          0          0

#just a list of samples that were actually used for phenotype associations (i.e. this is only used for subsetting)
df.clusters <- fread("code/finch_code/lowpass/15_GEMMA/file_input/Daphne_cluster1_samples.txt", header = F)

#this looks like this
#V1     V2
#1: 01Dap19024 fortis
#2: 01Dap19030 fortis
#3: 01Dap19032 fortis
#4: 01Dap19034 fortis
#5: 01Dap19035 fortis
#6: 01Dap19036 fortis

#bed region of the peaks
df.bed <- read.table(paste0("output/GEMMA/processed/", prefix, "_PEAKS.bed"))
names(df.bed) <- c("chr", "start","end", "peak")

#this looks like this
#chr    start      end peak
#1  chr1 61324213 61463374    1
#2  chr1 61954043 62619548    2
#3 chr1A 24867453 25022865    3
#4 chr1A 25267929 26113875    4
#5 chr1A 26603467 27142760    5
#6 chr1A 27804116 28390477    6

#remove positions if there is only one SNP
df.bed <- df.bed %>% 
  mutate(lg = end- start) %>% 
  filter(lg > 0)

#load the top p values (filtered > 0.1) to plot the signal with the other plots I generate below
df.logp <- read.table(paste0("output/GEMMA/processed/",prefix,"_OUTLIER_SNPs_TOP100.txt"))
#looks like:
#V1       V2       V3 V4
#1 chr1 61427453 7.792121  1
#2 chr1 61429631 7.454147  1
#3 chr1 61395560 7.158872  1
#4 chr1 61429326 7.112993  1
#5 chr1 61383660 7.075170  1
#6 chr1 61463374 7.061348  1

df.genes <- read.table(paste0("output/GEMMA/processed/",prefix,"_lmm_GENES.txt"), header = T)
#chr peak    start      end length num.snps mean.log_p max.log_p chr_ordered peak.original peak.new            gene_id strand gene_name
#1 chr1    1 61324213 61463374 139161       16   6.977724  7.792121        chr1             1        1 ENSCPVG00005007675      +       RB1
#2 chr1    1 61324213 61463374 139161       16   6.977724  7.792121        chr1             1        1 ENSCPVG00005007776      -     LPAR6
#3 chr1    1 61324213 61463374 139161       16   6.977724  7.792121        chr1             1        1 ENSCPVG00005007825      -    Rcbtb2
#4 chr1    2 61954043 62619548 665505      845   8.493228 13.100340        chr1             2        2 ENSCPVG00005008201      -    SPRYD7
#5 chr1    2 61954043 62619548 665505      845   8.493228 13.100340        chr1             2        2 ENSCPVG00005008213      +    TRIM13
#6 chr1    2 61954043 62619548 665505      845   8.493228 13.100340        chr1             2        2 ENSCPVG00005008220      +     kcnrg


###################################

#path must exist:
#output/GEMMA/genotype_plots/

#check to see if any NAs. Can cause problems later
min(colSums(is.na(df.genos[,-(1:4)])))
#could use this if there are
#df.genos.naomit <- na.omit(df.genos)

#empty lists to save plots and dataframes to
list.plots1 <- list() 
mega.list <- list() 

#three colors for three PCAs and NA (grey)
colors <- c("#440154FF","#2D708EFF","#FDE725FF", "grey")

#iterate for each peak in bed files

for (peak.idx in unique(df.bed$peak)){
  
  #take bed file, pull out peak of interest then filter by it
  df.bed.idx <- df.bed  %>% filter(peak == peak.idx)
  chr.idx <- unique(df.bed.idx$chr)
  print(paste0("running ", peak.idx, " on ", chr.idx))
  
  #subset geno df by the peak region 
  df.genos.sub <- df.genos %>% filter(CHROM == chr.idx & POS > min(df.bed.idx$start) & POS < max(df.bed.idx$end))
  
  if (nrow(df.genos.sub)> 0){
    
    # set up genotypes for plotting ------------------------------------------
    #only use individuals from the cluster of interest and remove position/chr colums

    df.genos.raw <- df.genos.sub %>% select(-CHROM, -REF, -ALT)
    #fread reads in columns with "XX" as charecters
    df.genos.raw.num <- sapply(df.genos.raw, as.numeric )
    
    # pca for groupings -------------------------------------------------------
    #if you only want to run PCA with the focal cluster ran with GWAS
    #df.x <- df.genos.sub %>% select(POS, df.clusters$V1) 
    
    df.x <- df.genos.raw
    
    row.names(df.x) <- df.x$POS ; df.x$POS <- NULL
    df.xt <- t(df.x)
    pca2 <- prcomp(df.xt, center = T, scale = T)
    
    df.pca2.values <- as.data.frame(pca2$x[,1:10])
    df.pca2 <- as.data.frame(summary(pca2)$importance) 
    
    #plot(df.pca2.values$PC2 ~ df.pca2.values$PC1)
    
    ##Next, run a k means cluster algorithm, using 25 samples as random, iterating up to a maximum of 100 times. 
    ##In these haplotypes, usually PC1 is enough to separate, but can consider if more is better
    k <- kmeans(df.pca2.values[,c(1)], 3, nstart = 25, iter.max = 100)
    
    ##In order to set a threshold of confidence use silhouette from the cluster package
    ##reasonable summary of this method here: https://towardsdatascience.com/clustering-analysis-in-r-using-k-means-73eca4fb7967
    ##values cluser to 0 are ambiguous, negative values are incorrect, positive values are better fit
    sil <- silhouette(k$cluster, dist(df.pca2.values[,c(1)]))
    
    ##Assign the cluster value and the confidence value to the main dataframe
    df.pca2.values$cluster <- paste("cluster",k$cluster, sep = "_")
    df.pca2.values$sil <- sil[,3]
    
    ##assign a threshold. In this case, values < 0.5 clearly cluster intermediately between homozygous and heterozygous individuals
    df.pca2.values <- df.pca2.values %>% mutate(cluster.pca = ifelse(sil < 0.5, NA, cluster))
    df.pca2.values$meaningful.unique <- row.names(df.pca2.values)
    
    df.meta <- V2.mdb[,c("meaningful.unique", "Species")]
    tmp <- left_join(df.pca2.values, df.meta, by = "meaningful.unique")
    ggplot(tmp)  + 
      geom_point(aes(x = PC1, y = PC2, color = Species))
    ggsave(paste0("output/GEMMA/genotype_plots/full_dataset_first_PCA_",peak.idx, ".pdf"), width = 8, height = 6)
    
    ggplot(tmp)  + 
      geom_point(aes(x = PC1, y = PC2, color = cluster.pca))
    ggsave(paste0("output/GEMMA/genotype_plots/full_dataset_first_color_cluster_PCA_",peak.idx, ".pdf"), width = 8, height = 6)
    
    
    df.geno.clusters.full <- left_join(V2.mdb, df.pca2.values[,c("meaningful.unique", "cluster.pca", "PC1")], by = "meaningful.unique")
    
    # assign clusters, this is just column naming ------------------------------------------------
    df.geno.clusters.full$cluster <- df.geno.clusters.full$cluster.pca ; df.geno.clusters.full$cluster.pca <- NULL
    #####THEN MAKE A DATAFRAME FOR ONLY PLOTTING THE CLUSTER OF INTEREST
    df.geno.clusters <- df.geno.clusters.full %>% filter(meaningful.unique %in% df.clusters$V1)
    ###########################################################################
    
    #get bill size mean for each cluster
    clust.mean <- df.geno.clusters %>% filter(!is.na(cluster)) %>% 
      group_by(cluster) %>% 
      summarise(bill.PC1 = mean(bill.PC1, na.rm = T),
                mean.PC1 = mean(PC1, na.rm = T)) %>% 
      arrange(mean.PC1) #%>% mutate(cluster = fct_reorder(clust.mean$cluster, clust.mean$bill.size))
    
    #ID heterozygote as the intermediate genotype in the PCA plot
    heterozygote.call <- clust.mean[2,1]$cluster
    
    #ID large and small homozygous alternative haplotypes (assuming relationship with beak size)
    df.large <- clust.mean %>% filter(cluster!=heterozygote.call) %>% arrange(-bill.PC1) %>% head(n=1) %>% select(cluster)
    homozygous.large.call <- df.large$cluster
    df.small <- clust.mean %>% filter(cluster!=heterozygote.call) %>% arrange(-bill.PC1) %>% tail(n=1) %>% select(cluster)
    homozygous.small.call <- df.small$cluster
    
    #I use clust mean later to arrange box plots so just re order it here
    clust.mean$cluster <- factor(clust.mean$cluster, levels = c(as.character(homozygous.small.call), as.character(heterozygote.call), as.character(homozygous.large.call)))
    clust.mean <- clust.mean %>% arrange(cluster)
    clust.mean$genotype <- c("AA","AB","BB")
    clust.mean$order <- 1:3
    clust.mean$cluster.label <- paste0(clust.mean$cluster, " (", clust.mean$genotype, ")")
    clust.mean$cluster.label <- fct_reorder(clust.mean$cluster.label, clust.mean$order)
    
    clust.match <- clust.mean[c("cluster", "genotype", "cluster.label")]
    
    df.geno.clusters <- left_join(df.geno.clusters, clust.match, by = "cluster") %>% 
      mutate(cluster = factor(cluster, levels = clust.match$cluster))
    df.geno.clusters.full <- left_join(df.geno.clusters.full, clust.match, by = "cluster") %>% 
      mutate(cluster = factor(cluster, levels = clust.match$cluster))
    
    # ID most differentiated SNPs ---------------------------------------------
    
    name.geno <- df.geno.clusters %>% select(meaningful.unique, genotype)
    
    df.genos.sub.num <- cbind(df.genos.sub[,c(1:4)], df.genos.raw.num)
    
    long.df.genos.sub <- df.genos.sub.num %>% 
      pivot_longer(cols = -c("CHROM", "POS", "ALT", "REF"), names_to = "sample", values_to = "genos") %>% 
      mutate(genos = ifelse(genos == "./." | is.na(genos), "NA", genos))
    
    long.df.genos.sub <- long.df.genos.sub %>% 
      left_join(name.geno, by = c('sample' = 'meaningful.unique'))
    
    topDAF <- long.df.genos.sub %>%
      filter(sample %in% df.clusters$V1) %>% 
      filter(genotype == "AA" | genotype == "BB") %>% 
      group_by(POS, genotype) %>% 
      summarise(sum.geno = sum(genos, na.rm = T), .groups = "drop") %>% 
      pivot_wider(values_from = sum.geno, names_from = genotype) %>% 
      mutate(AAF = AA / (AA+BB),
             BAF = BB / (AA+BB),
             DAF = abs(BAF - AAF)) %>% 
      ungroup() %>% 
      slice_max(DAF, n = 5)
    
    df.x <- df.genos.raw
    
    #only do top 5 DAF snps between alt haplotypes
    df.x <- df.x %>% filter(POS %in% topDAF$POS)
    
    row.names(df.x) <- df.x$POS ; df.x$POS <- NULL
    df.xt <- t(df.x)
    pca2 <- prcomp(df.xt, center = T, scale = T)
    
    df.pca2.values <- as.data.frame(pca2$x[,1:4])
    
    df.pca2 <- as.data.frame(summary(pca2)$importance) 
    
    #plot(df.pca2.values$PC2 ~ df.pca2.values$PC1)
    ##Next, run a k means cluster algorithm, using 25 samples as random, iterating up to a maximum of 100 times. 
    ##In these haplotypes, usually PC1 is enough to separate, but can consider if more is better
    k <- kmeans(df.pca2.values[,c(1)], 3, nstart = 25, iter.max = 100)
    
    ##In order to set a threshold of confidence use silhouette from the cluster package
    ##reasonable summary of this method here: https://towardsdatascience.com/clustering-analysis-in-r-using-k-means-73eca4fb7967
    ##values cluser to 0 are ambiguous, negative values are incorrect, positive values are better fit
    sil <- silhouette(k$cluster, dist(df.pca2.values[,c(1)]))
    
    ##Assign the cluster value and the confidence value to the main dataframe
    df.pca2.values$cluster <- paste("cluster",k$cluster, sep = "_")
    df.pca2.values$sil <- sil[,3]
    
    ##assign a threshold. In this case, values < 0.5 clearly cluster intermediately between homozygous and heterozygous individuals
    df.pca2.values <- df.pca2.values %>% mutate(cluster.pca = ifelse(sil < 0.5, NA, cluster))
    df.pca2.values$meaningful.unique <- row.names(df.pca2.values)
    
    
    # re run the seting of the AA and BB genos --------------------------------
    
    df.geno.clusters.full <- left_join(V2.mdb, df.pca2.values[,c("meaningful.unique", "cluster.pca", "PC1")], by = "meaningful.unique")
    
    # USE PCA OR DENDRO METHOD ------------------------------------------------
    ####TURN ON OR OFF USING PCA AS OUTPUT
    df.geno.clusters.full$cluster <- df.geno.clusters.full$cluster.pca ; df.geno.clusters.full$cluster.pca <- NULL
    #####THEN MAKE A DATAFRAME FOR ONLY PLOTTING THE CLUSTER OF INTEREST
    df.geno.clusters <- df.geno.clusters.full %>% filter(meaningful.unique %in% df.clusters$V1)
    ###########################################################################
    
    #get bill size mean for each cluster
    clust.mean <- df.geno.clusters %>% filter(!is.na(cluster)) %>% 
      group_by(cluster) %>% 
      summarise(bill.PC1 = mean(bill.PC1, na.rm = T),
                mean.PC1 = mean(PC1, na.rm = T)) %>% 
      arrange(mean.PC1) #%>% mutate(cluster = fct_reorder(clust.mean$cluster, clust.mean$bill.size))
    
    #ID heterozygote as the intermediate genotype in the PCA plot
    heterozygote.call <- clust.mean[2,1]$cluster
    
    #ID large and small homozygous alternative haplotypes (assuming relationship with beak size)
    df.large <- clust.mean %>% filter(cluster!=heterozygote.call) %>% arrange(-bill.PC1) %>% head(n=1) %>% select(cluster)
    homozygous.large.call <- df.large$cluster
    df.small <- clust.mean %>% filter(cluster!=heterozygote.call) %>% arrange(-bill.PC1) %>% tail(n=1) %>% select(cluster)
    homozygous.small.call <- df.small$cluster
    
    #I use clust mean later to arrange box plots so just re order it here
    clust.mean$cluster <- factor(clust.mean$cluster, levels = c(as.character(homozygous.small.call), as.character(heterozygote.call), as.character(homozygous.large.call)))
    clust.mean <- clust.mean %>% arrange(cluster)
    clust.mean$genotype <- c("AA","AB","BB")
    clust.mean$order <- 1:3
    clust.mean$cluster.label <- paste0(clust.mean$cluster, " (", clust.mean$genotype, ")")
    clust.mean$cluster.label <- fct_reorder(clust.mean$cluster.label, clust.mean$order)
    
    clust.match <- clust.mean[c("cluster", "genotype", "cluster.label")]
    
    df.geno.clusters <- left_join(df.geno.clusters, clust.match, by = "cluster") %>% 
      mutate(cluster = factor(cluster, levels = clust.match$cluster))
    df.geno.clusters.full <- left_join(df.geno.clusters.full, clust.match, by = "cluster") %>% 
      mutate(cluster = factor(cluster, levels = clust.match$cluster))
    
    
    
    # -------------------------------------------------------------------------
    
    
    df.out.cluster <- df.geno.clusters.full[,c("meaningful.unique","genotype")] 
    names(df.out.cluster) <- c("meaningful.unique", paste("genotype",chr.idx, peak.idx, sep = "_"))
    mega.list[[peak.idx]] <- df.out.cluster
    
    
    
    # run a LM  and make boxplot---------------------------------------------------------------
    m <- summary(lm(unlist(df.geno.clusters[,"bill.PC1"]) ~ unlist(df.geno.clusters[,"cluster"])))
    i<-chr.idx
    df.lm1 <- data.frame(NULL)
    df.lm1[i, "region"] <- "cluster"
    df.lm1[i, "r2"] <- round(m$adj.r.squared, 3)
    df.lm1[i, "pval"] <- m$coefficients[2,4] #pval
    df.lm1[i, "coef"] <- m$coefficients[2,1] #coef
    df.lm1[i, "n"] <- nrow(df.geno.clusters)
    
    df.geno.clusters.noNA <- df.geno.clusters %>% filter(!is.na(cluster))
    
    df.geno.clusters.noNA <- df.geno.clusters.noNA %>% filter(meaningful.unique %in% df.clusters$V1)
    
    list.plots1[[i]] <- ggplot() + geom_jitter(data = df.geno.clusters.noNA,  aes_string(x = "cluster.label", y = "bill.PC1"),color = "grey", width=0.15) +
      geom_boxplot(data = df.geno.clusters.noNA, aes_string(x = "cluster.label", y = "bill.PC1", fill = "cluster.label"),alpha = 0.5, outlier.shape = NA) + labs(subtitle = paste(i, peak.idx, sep = ":")) + theme_bw() +
      geom_text(aes(label = paste0("R2 = ",df.lm1$r2, ",p=",round(df.lm1$pval,4)), x = 2-0.5, y = max(df.geno.clusters.noNA$bill.PC1, na.rm = T) - 1)) +
      scale_fill_manual(values = colors) +
      theme(legend.position = "none") +
      labs(x = "Genotype", y = "Bill size (PC1)")
    
    p3 <- list.plots1[[i]] 
    #ggsave(paste0("output/GEMMA/Daphne_cluster1_boxplot_",chr,".png"), width = 8, height = 8)
    
    m <- summary(lm(unlist(df.geno.clusters[,"weight"]) ~ unlist(df.geno.clusters[,"cluster"])))
    i<-chr.idx
    df.lm3 <- data.frame(NULL)
    df.lm3[i, "region"] <- "cluster"
    df.lm3[i, "r2"] <- round(m$adj.r.squared, 3)
    df.lm3[i, "pval"] <- m$coefficients[2,4] #pval
    df.lm3[i, "coef"] <- m$coefficients[2,1] #coef
    df.lm3[i, "n"] <- nrow(df.geno.clusters)
    
    df.geno.clusters.noNA <- df.geno.clusters %>% filter(!is.na(cluster))
    
    df.geno.clusters.noNA <- df.geno.clusters.noNA %>% filter(meaningful.unique %in% df.clusters$V1)
    
    p.weight <- ggplot() + geom_jitter(data = df.geno.clusters.noNA,  aes_string(x = "cluster.label", y = "weight"),color = "grey", width=0.15) +
      geom_boxplot(data = df.geno.clusters.noNA, aes_string(x = "cluster.label", y = "weight", fill = "cluster.label"),alpha = 0.5, outlier.shape = NA) + labs(subtitle = paste(i, peak.idx, sep = ":")) + theme_bw() +
      geom_text(aes(label = paste0("R2 = ",df.lm3$r2), x = 2-0.5, y = max(df.geno.clusters.noNA$weight, na.rm = T) - 1)) +
      scale_fill_manual(values = colors) +
      theme(legend.position = "none") +
      xlab("Genotype")
    
    
    # run a LM  and make boxplot - beak shape ---------------------------------------------------------------
    df.geno.clusters$bill.shape <- df.geno.clusters$bill.length / df.geno.clusters$bill.depth
    m <- summary(lm(unlist(df.geno.clusters[,"bill.PC2"]) ~ unlist(df.geno.clusters[,"cluster"])))
    i<-chr.idx
    df.lm2 <- data.frame(NULL)
    df.lm2[i, "region"] <- "cluster"
    df.lm2[i, "r2"] <- round(m$adj.r.squared, 3)
    df.lm2[i, "pval"] <- m$coefficients[2,4] #pval
    df.lm2[i, "coef"] <- m$coefficients[2,1] #coef
    df.lm2[i, "n"] <- nrow(df.geno.clusters)
    
    df.geno.clusters.noNA <- df.geno.clusters %>% filter(!is.na(cluster))
    p7  <- ggplot() + geom_jitter(data = df.geno.clusters.noNA,  aes_string(x = "cluster.label", y = "bill.PC2"),color = "grey", width=0.15) +
      geom_boxplot(data = df.geno.clusters.noNA, aes_string(x = "cluster.label", y = "bill.PC2", fill = "cluster.label"), outlier.shape = NA, alpha = 0.5) + labs(subtitle = paste(i, peak.idx, sep = ":")) + theme_bw() +
      geom_text(aes(label = paste0("R2 = ",df.lm2$r2, ",p=",round(df.lm2$pval,4)), x = 2-0.5, y = max(df.geno.clusters.noNA$bill.PC2, na.rm = T) - .25)) +
      scale_fill_manual(values = colors) +
      theme(legend.position = "none") +
      labs(x = "Genotype", y = "Bill shape (PC1)")

    p3 + p7
    ggsave(paste0("output/GEMMA/genotype_plots/bill_size_shape_boxplot_",chr.idx, "_",peak.idx, ".pdf"), width = 12, height = 6)
    
    # plot pca cluster ------------------------------------------------------------
    df.pca2.values <-  df.pca2.values %>% 
      mutate(cluster = factor(cluster.pca, levels = clust.match$cluster)) 
    df.pca2.values <- left_join(df.pca2.values, clust.match, by = "cluster")
    
    p.pca <- df.pca2.values %>% 
      filter(meaningful.unique %in% df.clusters$V1) %>% 
      ggplot(aes(x = PC1, y = PC2, fill = cluster.label), color = "black") + 
      geom_point(pch = 21, size = 3) + 
      xlab(paste0("PC1 (", 100 * df.pca2[2,1],"%)")) +
      ylab(paste0("PC2 (", 100 * df.pca2[2,2],"%)")) +
      theme_bw() + scale_fill_manual(values = colors)
    
    
    df.pca2.values <- left_join(df.pca2.values, V2.mdb[c("meaningful.unique", "Species")], by = "meaningful.unique")
    
    
    # pca histogram -----------------------------------------------------------
    p.histpca <- df.pca2.values %>% 
      filter(meaningful.unique %in% df.clusters$V1) %>% 
      ggplot() + 
      geom_histogram(aes(x = PC1, fill = cluster.label), bins = 30) + 
      theme_bw() + scale_fill_manual(values = colors)
    
    # plot PCA full dataset ---------------------------------------------------
    p.pca.full <- df.pca2.values %>% ggplot(aes(x = PC1, y = PC2, fill = cluster.label), color = "black") + 
      geom_point(pch = 21, size = 3) + 
      xlab(paste0("PC1 (", 100 * df.pca2[2,1],"%)")) +
      ylab(paste0("PC2 (", 100 * df.pca2[2,2],"%)")) +
      theme_bw() + scale_fill_manual(values = colors)
    p.pca.full
    ggsave(paste0("output/GEMMA/genotype_plots/full_dataset_PCA_",peak.idx, ".pdf"), width = 8, height = 6)
    
    # heatmap -----------------------------------------------------------------
    
    df.cluster.map <- df.geno.clusters %>% select(meaningful.unique, cluster, cluster.label)
  
    df.genos.sub.num <- cbind(df.genos.sub[,c(1:4)], df.genos.raw.num)
    
    long.df.genos.sub <- df.genos.sub.num %>% 
      pivot_longer(cols = -c("CHROM", "POS", "ALT", "REF"), names_to = "sample", values_to = "genos") %>% 
      mutate(genos = ifelse(genos == "./." | is.na(genos), "NA", genos))
    #colors <- c('pink','black','yellow', "white")
    
    #only do those 5 SNPs used
    long.df.genos.sub <- long.df.genos.sub %>% filter(POS %in% topDAF$POS)
    
    long.df.genos.sub <- long.df.genos.sub %>% 
      filter(sample %in% df.clusters$V1) %>% #cluster map is only the fortis samples
      left_join(df.cluster.map, by = c('sample' = 'meaningful.unique'))
    
    region.lg <- (max(long.df.genos.sub$POS) - min(long.df.genos.sub$POS)) / 1000
    
    #polarize SNPs. First get one individual from the first cluster
    ref.ind <- long.df.genos.sub %>% filter(cluster == clust.mean[1,1]$cluster) %>% distinct(sample) %>% head(n = 1)
    long.df.genos.sub <- long.df.genos.sub %>% filter(sample == ref.ind$sample) %>% 
      select(POS, genos, ref.geno = genos) %>% 
      right_join(long.df.genos.sub, by = "POS") %>%
      mutate(geno.polarized = ifelse(genos == 0 & ref.geno == 0, 0, #if ref is 0, alt is 0
                                     ifelse(genos == 2 & ref.geno == 2, 0, #if ref is alt, alt is 0
                                            ifelse(genos == 0 & ref.geno == 2, 2,# if ref is alt, alt is 2
                                                   ifelse(genos == 2 & ref.geno == 0, 2, #if ref is 0, alt is 2
                                                          ifelse(genos == 1, 1, #leave hets alone
                                                                 ifelse(ref.geno == 1, genos, NA))))))) #if het just keep OG geno
    
    p4 <- long.df.genos.sub %>% 
      #filter(!is.na(cluster)) %>% 
      ggplot(aes(x = as.factor(POS), y = sample)) + 
      geom_tile(aes(fill = as.factor(geno.polarized))) +
      scale_fill_manual(values=c("#fc8d59", "black","#91bfdb"))+
      #theme(axis.text.x = element_text(angle = 90)) +
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            legend.position = "none") +
      facet_wrap(~cluster.label, scale="free_y", nrow=4) +
      ggtitle(paste0("region size (kb) = ", region.lg)) +
      xlab(NULL)
    
    ##cluster frequency by year -----------------------------------------------
    ##based on year born
    df.cluster.all <- df.pca2.values %>% select(meaningful.unique, cluster.label)
    
    df.cluster.map1 <- df.cluster.all %>% left_join(V2.mdb[,c("meaningful.unique", "Species","First.year.min", "Last.year")], by = "meaningful.unique")
    
    # cluster frequency by year -----------------------------------------------  
    #based on year alive
    
    df.year.spread <- df.cluster.map1 %>% 
      filter(!is.na(First.year.min)) %>% filter(!is.na(Last.year)) %>% 
      select(meaningful.unique, cluster.label, First.year.min, Last.year, Species) %>% 
      mutate(year = map2(First.year.min, Last.year, `:`)) %>% 
      select(-First.year.min, -Last.year) %>% 
      unnest(cols = c(year))
    
    df.sum.year2 <- df.year.spread %>% group_by(year, Species) %>% 
      summarise(year.n = n(), .groups = "keep") 
    
    df.year.prop2 <- df.year.spread %>% 
      group_by(year, cluster.label, Species) %>% 
      summarise(geno_sum = n(), .groups = "keep") %>% 
      left_join(df.sum.year2, by = c("year", "Species")) %>% 
      mutate(prop = geno_sum / year.n)
    
    p6 <- df.year.prop2 %>% filter(year > 1987 & Species == "fortis") %>% 
      filter(!is.na(cluster.label)) %>% 
      ggplot(aes(x = year, y = prop, color = cluster.label)) +
      geom_point(aes(size = year.n)) + geom_line() +
      scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1986, 2012, by = 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_blank(),
            legend.position = "bottom") +
      scale_color_manual(values = colors) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ylab("Frequency of genotype")
    
    p6
    ggsave(paste0("output/GEMMA/genotype_plots/fortis_year_change_",peak.idx, ".pdf"), width = 8, height = 6)
    
    p8 <- df.year.prop2 %>% filter(year > 1987 & Species == "scandens") %>% 
      filter(!is.na(cluster.label)) %>% 
      ggplot(aes(x = year, y = prop, color = cluster.label)) +
      geom_point(aes(size = year.n)) + geom_line() +
      scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1986, 2012, by = 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_blank(),
            legend.position = "bottom") +
      scale_color_manual(values = colors) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ylab("Frequency of genotype")
    
    p8
    ggsave(paste0("output/GEMMA/genotype_plots/scandens_year_change_",peak.idx, ".pdf"), width = 8, height = 6)
    
    
    # plot all together -------------------------------------------------------
    
    #(p1 + ~ plot(df.genos.dend)) / (p3 + p4) / p5
    (p.histpca + p4) / (p3 + p7) / p6 / p8
    
    ggsave(paste0("output/GEMMA/genotype_plots/Daphne_cluster1_multi",chr.idx,"_",peak.idx, ".png"), width = 8.5, height = 17)
    
    
    
  }
}

mega.list[sapply(mega.list, is.null)] <- NULL

df.finch.clusters <- Reduce(function(...) merge(..., by=c("meaningful.unique"), all.x=TRUE), mega.list)

write_csv(df.finch.clusters,"output/GEMMA/Daphne_cluster1_all_haplotypes.csv")

df.metadata <- V2.mdb %>% select(meaningful.unique,TubeNum, AltNum.08, Page.01, Slot.02, Species, Island, Symbol, contains("Year"), starts_with("bill"), weight, tarsus, wing, N.broods, 
                                        nufledg, Sex, Recruits.Total, Recruits.male, Recruits.female, Longevity.yrs, mother, EP.father, social.father, Year.sampled)

df.finch.clusters <- df.finch.clusters %>% left_join(df.metadata, by = "meaningful.unique")
df.finch.clusters <- df.finch.clusters %>% mutate(bill.size = bill.length + bill.width + bill.depth,
                                          bill.shape = bill.length / bill.depth)


#identify corresponding regions from original analysis -----------------

#confirmed these are the same strong associations by examing the genotype analysis

df.bed.og <- read.table("output/GEMMA/processed/autosomes_Daphne_cluster1_multivariate_INTs_lmm1_PEAKS.bed")

six.peaks <- c(1,9,16,19,21,23)

#9 is hmga2, 16 is alx1
df.bed6 <- df.bed.og %>% 
  filter(V4 %in% six.peaks)


library(GenomicRanges)
a <- GRanges(df.bed6$V1, IRanges(as.numeric(df.bed6$V2), as.numeric(df.bed6$V3)))
b <- GRanges(df.bed$chr, IRanges(as.numeric(df.bed$start), as.numeric(df.bed$end)))

comp_overlap <- findOverlaps(b, a)
comp.annotated <- cbind(df.bed[queryHits(comp_overlap), ], df.bed6[subjectHits(comp_overlap), ])

strong_assoc <- paste("genotype", comp.annotated$chr, comp.annotated$peak, sep = "_")

# -------------------------------------------------------------------------

other_cols <- grep("genotype_", names(df.finch.clusters), value = T, invert = T)

df.finch.clusters.strong <- df.finch.clusters %>% 
  select(other_cols, strong_assoc)

write_csv(df.finch.clusters.strong, "output/GEMMA/All_lowpass_samples_strong_associations.csv", na = "")

df.finch.clusters.weak <- df.finch.clusters %>% 
  select(-strong_assoc)

write_csv(df.finch.clusters.weak, "output/GEMMA/All_lowpass_samples_weak_associations.csv", na = "")




















# correlations ------------------------------------------------------------
#for cramers v
library(DescTools)
#for ggplot correlation
library(ggcorrplot)
library(gtools) # combination
library(data.table) # data mgmt

tmp.df <- df.finch.clusters %>% select(cluster_chr1_10, cluster_chr1A_53, cluster_chr1A_55, cluster_chr2_133, cluster_chr9_381)
#CramerV(table(tmp.df))

model.matrix(~0+., data=tmp.df) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)
ggsave("output/GEMMA/genotype_plots/correl_genotype_clusters.png", width = 6, height = 6)


#cramers V from: https://github.com/Julien-Yacine/cramer_v_heatmap
cv.test = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
              (length(x)[1] * (min(length(unique(x))[1],length(unique(y))[1]) - 1)))
  return(as.numeric(CV))
}

# Apply the function to the combination of categorical variable
v_cramer_all <- function(cat_var, df){
  cat_var_grid <- data.table(combinations(n = length(cat_var), r = 2, v = cat_var, repeats.allowed = F))
  
  do.call(rbind,
          apply(cat_var_grid, 1, function(x){
            tmp <- as.character(x)
            vec1 <- base::unlist(df[,tmp[1]])
            vec2 <- base::unlist(df[,tmp[2]])
            
            data.table(
              variable_x = tmp[1],
              variable_y = tmp[2],
              chi2 = chisq.test(x = vec1, vec2, correct=FALSE)$p.value,
              v_cramer = cv.test(x = vec1, y = vec2)
            )
          }))
  
}

results <- v_cramer_all(cat_var = colnames(tmp.df), df = tmp.df)
p.cramers <- ggplot(results, aes(variable_x, variable_y)) +
  geom_tile(aes(fill = v_cramer), colour = "black") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("Cramer's V heatmap")
p.cramers
ggsave("output/GEMMA/genotype_plots/cramersV_genotype_clusters.png", width = 6, height = 6)
