library(noia)
library(tidyverse)
library(effectsize)
library(patchwork)
#for pretty legends:
library(cowplot)
library(grid)
library(gridExtra)
library(ggtext)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

#setup data. 
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

#drop outliers (two fortis immigrants, one may be maginrostris)
df.all.genos <- df.all.genos %>% 
  filter(meaningful.unique!="01Dap21269" & meaningful.unique!="01Dap21274")


#update sex with field sex if genetic is missing
df.all.genos <- df.all.genos %>% 
  mutate(genetic.sex = as.factor(gsub("unassigned", NA, genetic.sex)),
         Sex2 = ifelse(Sex == 1, "Male", 
                       ifelse(Sex ==2, "Female", NA))) %>% 
  mutate(munch.sex = ifelse(is.na(genetic.sex) & !is.na(Sex2), as.character(Sex2), as.character(genetic.sex)))

df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")

df.dap.genos.num <- df.dap.genos 
df.dap.genos.num[df.dap.genos.num=="AB"]<-2
df.dap.genos.num[df.dap.genos.num=="AA"]<-1
df.dap.genos.num[df.dap.genos.num=="BB"]<-3

gwas_loci <- c("gwas_genotype_chr1_2_ext", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18_ext", 
               "gwas_genotype_chr9_20")

asm28loci <- grep("asm28loci", names(df.all.genos), value = T)

asm28loci.minus.linked <- grep(paste(c("_2$","_4","_5","_6"), collapse="|"), asm28loci, invert=T, value=T)
asm28loci.minus.gwas <- grep(paste(c("_1$","_3","_7","_27"), collapse="|"), asm28loci.minus.linked, invert=T, value=T)

all_loci <- c(gwas_loci, asm28loci.minus.gwas)

#hand curated by LD profile and highest pvalue SNP in fortis gwas
ld_filt_loc <- c(gwas_loci, "asm28loci_chr2_10", "asm28loci_chr2_12","asm28loci_chr2_17",
                 "asm28loci_chr2_19","asm28loci_chr2_20","asm28loci_chr3_22", "asm28loci_chr3_24", "asm28loci_chr5_25",
                 "asm28loci_chr7_26", "asm28loci_chr25_28")

##tibble(locus = ld_filt_loc) %>% write_csv("output/GEMMA/effect_sizes/ld_prune_locus_names.csv")

ld_filt_loc.names <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")


# include ancestry as predictor -------------------------------------------

calc_es_with_ancestry <- function(species, target_loci, output.name, cluster = NA, anc_group1, anc_group2){
  
  df.finch <- df.dap.genos.num %>% filter(Species == species)
  
  if (is.na(cluster)){
    df.finch <- df.dap.genos.num %>% filter(Species == species)  
  }else{
    df.finch <- df.dap.genos.num %>% filter(grm_cluster == cluster)
  }
  
  if(!is.na(cluster)){
    species <- cluster #for naming purposes
  }
  
  df.finch.sub <- df.finch %>% 
    select("weight","bill.length", "bill.depth", "bill.width",
           all_of(target_loci), "anc_fort", "anc_scan", "anc_mag", "anc_fuli")
  
  rownames(df.finch.sub) <- df.finch$meaningful.unique
  df.finch.sub <- df.finch.sub %>% 
    filter(!is.na(bill.width) & !is.na(bill.depth) & !is.na(bill.length))
  
  # calculate PCs only for species filtered  --------------------------------
  
  beak.pca <- prcomp(df.finch.sub[,c("bill.depth", "bill.width", "bill.length")], center = T, scale = T)
  
  df.finch.sub <- cbind(df.finch.sub,beak.pca$x[,1:3])
  #replace existing values (which was built from 4 species)
  df.finch.sub$bill.PC1 <- df.finch.sub$PC1
  df.finch.sub$bill.PC2 <- df.finch.sub$PC2
  df.finch.sub$bill.PC3 <- df.finch.sub$PC3
  
  
  # orient PCs so that they are positive correlated with beak propor --------
  corr_PC1 <- cor(df.finch.sub$bill.PC1, 
                  df.finch.sub$bill.depth + df.finch.sub$bill.length + df.finch.sub$bill.width)
  
  corr_PC2 <- cor(df.finch.sub$bill.PC2, 
                  df.finch.sub$bill.length / df.finch.sub$bill.depth)
  
  df.finch.sub <- df.finch.sub %>% 
    mutate(bill.PC1 = ifelse(rep(corr_PC1, nrow(df.finch.sub)) < 0, bill.PC1*-1, bill.PC1),
           bill.PC2 = ifelse(rep(corr_PC2, nrow(df.finch.sub)) < 0, bill.PC2*-1, bill.PC2))
  
  # plot effect sizes -------------------------------------------------------
  
  #trick is that if not enough var models wont run
  #make list of vars with more than 2 
  #because of NAs included, needs to be more than 3
  keep_list <- sapply(lapply(df.finch.sub, unique), length) > 3
  
  df.tmp <- df.finch.sub[, keep_list]
  target_loci2 <- names(df.tmp)[names(df.tmp) %in% target_loci]
  
  testvars <- c(target_loci2, anc_group1, anc_group2)
  lmvars <- paste(testvars, collapse=" + ")

  m.PC1 <- lm(paste("bill.PC1 ~ ", lmvars), data = df.finch.sub, na.action = na.exclude)
  A.PC1 <- car::Anova(m.PC1, type = 3, singular.ok = T)
  eta.PC1 <- eta_squared(A.PC1, partial = F, alternative = "two.sided")
  
  #summary(lm(bill.PC1 ~ anc_scan + anc_fuli, data = df.finch.sub, na.action = na.exclude))
  
  m.PC2 <- lm(paste("bill.PC2 ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  A.PC2 <- car::Anova(m.PC2, type = 3, singular.ok = T)
  eta.PC2 <- eta_squared(A.PC2, partial = F, alternative = "two.sided")
  
  #plot(df.finch.sub$anc_scan, df.finch.sub$bill.PC2.sp)
  
  m.weight <- lm(paste("weight ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  A.weight <- car::Anova(m.weight, type = 3, singular.ok = T)
  eta.weight <- eta_squared(A.weight, partial = F, alternative = "two.sided")
  
  # add sign to es ----------------------------------------------------------
  #get coefficients for PC1 and PC2
  df.coef <- tibble(locus.mod = names(m.PC1$coefficients),
                    coef = m.PC1$coefficients,
                    coef2 = m.PC2$coefficients,
                    coef3 = m.weight$coefficients)
  
  df.coef$locus <- str_sub(df.coef$locus.mod,1,nchar(df.coef$locus.mod)-1)
  
  #fix ancestry naming
  df.coef$locus <- gsub("anc_for", "anc_fort", gsub("anc_sca", "anc_scan", gsub("anc_ful", "anc_fuli", df.coef$locus)))
  
  df.coef.sign <- df.coef %>% 
    filter(locus!="(Intercept") %>% 
    group_by(locus) %>% 
    summarise(sum.coef = sum(coef),
              sum.coef2 = sum(coef2),
              sum.coef3 = sum(coef3)) %>% 
    mutate(sign = ifelse(sum.coef >= 0, "+",
                         ifelse(sum.coef <0, "-", NA))) %>% 
    mutate(sign2 = ifelse(sum.coef2 >= 0, "+",
                          ifelse(sum.coef2 <0, "-", NA))) %>% 
    mutate(sign3 = ifelse(sum.coef3 >= 0, "+",
                          ifelse(sum.coef3 <0, "-", NA))) %>% 
    select(-sum.coef, -sum.coef2, -sum.coef3)
  # PCs ---------------------------------------------------------------------
  
  df.effect.sizePCs <- tibble(`Bill PC1` = eta.PC1$Eta2,
                              `Bill PC1 low` = eta.PC1$CI_low,
                              `Bill PC1 high` = eta.PC1$CI_high,
                              `Bill PC2` = eta.PC2$Eta2,
                              `Bill PC2 low` = eta.PC2$CI_low,
                              `Bill PC2 high` = eta.PC2$CI_high,
                              `Weight` = eta.weight$Eta2,
                              `Weight low` = eta.weight$CI_low,
                              `Weight high` = eta.weight$CI_high,
                              `Bill PC1 p` = A.PC1$`Pr(>F)`[2:(length(testvars)+1)],
                              `Bill PC2 p` = A.PC2$`Pr(>F)`[2:(length(testvars)+1)],
                              `Weight p` = A.weight$`Pr(>F)`[2:(length(testvars)+1)],
                              locus = testvars,
                              `Bill PC1 var` = rep(var(df.finch.sub$bill.PC1), length(testvars)),
                              `Bill PC2 var` = rep(var(df.finch.sub$bill.PC2), length(testvars))
  ) 

  df.effect.sizePCs <- df.effect.sizePCs%>% 
    pivot_longer(cols = -locus, values_to = "parameter", names_to = "phenotype") %>% 
    mutate(locus = factor(locus, levels = testvars))
  
  
  # add variation? ----------------------------------------------------------
  
  df.finch.sub2 <- df.finch.sub
  df.finch.sub2 <- df.finch.sub2 %>% 
    select(all_of(target_loci))
  df.finch.sub2$index <- 1:nrow(df.finch.sub2)
  
  #df.finch.sub2$genotype <- df.finch.sub2$HMGA2.simple
  
  df.finch.sub2[df.finch.sub2=="2"]<-"AB"
  df.finch.sub2[df.finch.sub2=="1"]<-"AA"
  df.finch.sub2[df.finch.sub2=="3"]<-"BB"
  
  df.finch.sub2.long <- df.finch.sub2 %>% 
    pivot_longer(cols = -index, values_to = "genotype", names_to = "locus") %>% 
    select(-index)
  
  
  df.freq <- df.finch.sub2.long %>% group_by(locus, genotype) %>% 
    summarise(count = n(), .groups = "keep")
  
  df.freq.haps <- df.freq %>% 
    pivot_wider(values_from = count, names_from = genotype) %>% 
    replace(is.na(.), 0) %>% 
    mutate(n = AA + AB + BB+ `NA`) %>% 
    mutate(af.A = ((2*AA) + AB)/(2*(n-`NA`)),
           af.B = ((2*BB) + AB)/(2*(n-`NA`)),
           q = AA / (n - `NA`), #genotype frequency of AA
           p = BB / (n - `NA`),
           af.sum = af.A + af.B,
           MAF = min(af.A, af.B),
           bAF = af.B,
           minGF = min(q, p))
  
   df.effect.sizePCs <- left_join(df.effect.sizePCs, df.freq.haps[,c("locus","MAF","bAF", "minGF")], by = "locus")
  
  miss_af <- df.freq.haps %>% filter(!locus %in% target_loci2) %>% select("locus","MAF","bAF", "minGF")
  # -------------------------------------------------------------------------
  
  
  #update sign iteratively
  df.effect.sizePCs <- left_join(df.effect.sizePCs, df.coef.sign, by = "locus")
  
  #update sign PC1
  df.effect.sizePCs_A <- df.effect.sizePCs %>% 
    filter(grepl("PC1", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign == "-", parameter * -1, parameter))
  
  #update sign PC2
  df.effect.sizePCs_B <- df.effect.sizePCs %>% 
    filter(grepl("PC2", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign2 == "-", parameter * -1, parameter))
  
  #update sign weight
  df.effect.sizePCs_C <- df.effect.sizePCs %>% 
    filter(grepl("Weight", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign3 == "-", parameter * -1, parameter))
  
  #extract p values
  df.effect.sizePCs_p <- df.effect.sizePCs %>% 
    filter(grepl(" p", phenotype) |  grepl("var", phenotype)) 
  
  #recombine 3 dfs
  df.effect.sizePCs <- rbind(df.effect.sizePCs_A, df.effect.sizePCs_B, df.effect.sizePCs_C, df.effect.sizePCs_p) %>% 
    select(-sign, -sign2, -sign3)
  
  max_val <- df.effect.sizePCs %>% filter
  
  p2 <- df.effect.sizePCs %>% 
    filter(!grepl("high", phenotype) & !grepl("low", phenotype)) %>% 
    filter(!grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    ggplot(aes(label = round(parameter, 2),x = locus, y = phenotype, fill = parameter)) + 
    geom_tile() +
    geom_text() +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  p2
  ggsave(paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_long_no_residuals_ancestry.pdf"), width = 11, height = 8.5)
  
  list.plots <- list()
  
  for (genotype in target_loci){
    
    list.plots[[genotype]] <- df.finch.sub %>% 
      filter(!is.na(get(genotype))) %>% 
      ggplot(aes_string(x = genotype, y = "bill.PC1")) +
      geom_jitter(width = 0.2) + 
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      theme_bw()
  }
  
  pZ <- wrap_plots(list.plots, nrow = 3)
  pZ
  ggsave(paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_", output.name,"_boxplot_no_residuals_ancestry.pdf"), width = 11, height = 8.5)
  
  df.effect.sizePCs.wide <- df.effect.sizePCs %>% 
    pivot_wider(values_from = parameter, names_from = phenotype) %>% 
    mutate(across(where(is.numeric), round, 3)) 
  
  #make null data for loci that were not variable enough to be included (set to 0 for pretty plotting)
  
  miss_loci <- target_loci[!target_loci %in% target_loci2]
  
  df_null <- tibble(MAF = miss_af$MAF,
                    bAF = miss_af$bAF,
                    minGF = miss_af$minGF,
                    `Bill PC1` = 0,
                    `Bill PC1 low` = 0,
                    `Bill PC1 high` = 0,
                    `Bill PC2` = 0,
                    `Bill PC2 low` = 0,
                    `Bill PC2 high` = 0,
                    `Weight` = 0,
                    `Weight low` = 0,
                    `Weight high` = 0,
                    `Bill PC1 p` = NA,
                    `Bill PC2 p` = NA,
                    `Weight p` = NA,
                    locus = miss_af$locus,
                    `Bill PC1 var` = rep(var(df.finch.sub$bill.PC1), nrow(miss_af)),
                    `Bill PC2 var` = rep(var(df.finch.sub$bill.PC2), nrow(miss_af))
  ) 
  
  df.effect.sizePCs.wide <- rbind(df.effect.sizePCs.wide, df_null)
  
  df.effect.sizePCs.wide$`Bill PC1 p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Bill PC1 p`, method = "BH")
  df.effect.sizePCs.wide$`Bill PC2 p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Bill PC2 p`, method = "BH")
  df.effect.sizePCs.wide$`Weight p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Weight p`, method = "BH")
  
  write_csv(df.effect.sizePCs.wide, paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_es_no_residuals_with_ancestry.csv"))
  
}

calc_es_with_ancestry("fortis", gwas_loci, "gwas_loci_with_ancestry", "cluster1", "anc_scan", "anc_fuli")
calc_es_with_ancestry("scandens", gwas_loci, "gwas_loci_with_ancestry", "cluster2", "anc_fort", "anc_fuli")

# -------------------------------------------------------------

plot_effect_size <- function(species, target_loci, output.name, cluster = NA, ancfilt = NA, resid = FALSE, sex = TRUE){
  
  df.finch <- df.dap.genos.num %>% filter(Species == species)
  
  if (is.na(cluster)){
    df.finch <- df.dap.genos.num %>% filter(Species == species)  
  }else{
    df.finch <- df.dap.genos.num %>% filter(grm_cluster == cluster)
  }
  
  if(!is.na(ancfilt) & species == "fortis"){
    #cutoff for ancestry is 25% (so top 75% quantile)
    cutoff <- quantile(df.finch$anc_fort, na.rm = T)[2]
    df.finch <- df.finch %>% 
      filter(anc_fort > cutoff)
    
    output.name <- paste(output.name, "ancfilt", sep = "_")
  }
  
  if(!is.na(ancfilt) & species == "scandens"){
    #cutoff for ancestry is 25% (so top 75% quantile)
    cutoff <- quantile(df.finch$anc_scan, na.rm = T)[2]
    df.finch <- df.finch %>% 
      filter(anc_scan > cutoff)
    output.name <- paste(output.name, "ancfilt", sep = "_")
  }
  
  if(!is.na(cluster)){
    species <- cluster #for naming purposes
  }
  
  df.finch.sub <- df.finch %>% 
    select("weight","bill.length", "munch.sex", "bill.depth", "bill.width",
           all_of(target_loci))
  
  rownames(df.finch.sub) <- df.finch$meaningful.unique
  df.finch.sub <- df.finch.sub %>% 
    filter(!is.na(bill.width) & !is.na(bill.depth) & !is.na(bill.length))
  
  # calculate PCs only for species filtered  --------------------------------
  
  beak.pca <- prcomp(df.finch.sub[,c("bill.depth", "bill.width", "bill.length")], center = T, scale = T)
  write_csv(data.frame(summary(beak.pca)$rotation), paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_importance.csv"))
  
  df.finch.sub <- cbind(df.finch.sub,beak.pca$x[,1:3])
  #replace existing values (which was built from 4 species)
  df.finch.sub$bill.PC1 <- df.finch.sub$PC1
  df.finch.sub$bill.PC2 <- df.finch.sub$PC2
  df.finch.sub$bill.PC3 <- df.finch.sub$PC3
  

  # orient PCs so that they are positive correlated with beak propor --------

  
  corr_PC1 <- cor(df.finch.sub$bill.PC1, 
            df.finch.sub$bill.depth + df.finch.sub$bill.length + df.finch.sub$bill.width)
  
  corr_PC2 <- cor(df.finch.sub$bill.PC2, 
                        df.finch.sub$bill.length / df.finch.sub$bill.depth)
  
  df.finch.sub <- df.finch.sub %>% 
    mutate(bill.PC1 = ifelse(rep(corr_PC1, nrow(df.finch.sub)) < 0, bill.PC1*-1, bill.PC1),
           bill.PC2 = ifelse(rep(corr_PC2, nrow(df.finch.sub)) < 0, bill.PC2*-1, bill.PC2))
  
  df.finch.sub %>% 
    ggplot() + geom_point(aes(x = bill.PC2, y = bill.length/bill.depth))
  
  ggsave(paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_PC2_len_dep.pdf"), 
         width = 8, height = 6)
  
  if(resid == TRUE & sex == TRUE){
    
    df.finch.sub$bill.PC1 <- resid(lm(bill.PC1 ~ weight + munch.sex, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.PC2 <- resid(lm(bill.PC2 ~ weight + munch.sex, data = df.finch.sub, na.action = na.exclude))
    
    df.finch.sub$bill.length <- resid(lm(bill.length ~ weight + munch.sex, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.depth <- resid(lm(bill.depth ~ weight + munch.sex, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.width <- resid(lm(bill.width ~ weight + munch.sex, data = df.finch.sub, na.action = na.exclude))
  }
  
  if(resid == TRUE & sex == FALSE){
    
    df.finch.sub$bill.PC1 <- resid(lm(bill.PC1 ~ weight, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.PC2 <- resid(lm(bill.PC2 ~ weight, data = df.finch.sub, na.action = na.exclude))
    
    df.finch.sub$bill.length <- resid(lm(bill.length ~ weight, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.depth <- resid(lm(bill.depth ~ weight, data = df.finch.sub, na.action = na.exclude))
    df.finch.sub$bill.width <- resid(lm(bill.width ~ weight, data = df.finch.sub, na.action = na.exclude))
  }
  
  # plot effect sizes -------------------------------------------------------
  
  #trick is that if not enough var models wont run
  #make list of vars with more than 2 
  #because of NAs included, needs to be more than 3
  keep_list <- sapply(lapply(df.finch.sub, unique), length) > 3

  df.tmp <- df.finch.sub[, keep_list]
  target_loci2 <- names(df.tmp)[names(df.tmp) %in% target_loci]
  
  lmvars <- paste(target_loci2, collapse=" + ")
  
  m.depth <- lm(paste("bill.depth ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)

  eta.depth <- eta_squared(car::Anova(m.depth, type = 3, singular.ok = T), partial = F)
  
  m.width <- lm(paste("bill.width ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  
  eta.width <- eta_squared(car::Anova(m.width, type = 3, singular.ok = T), partial = F)
  
  m.length <- lm(paste("bill.length ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  
  eta.length <- eta_squared(car::Anova(m.length, type = 3,singular.ok = T), partial = F)
  
  m.PC1 <- lm(paste("bill.PC1 ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  A.PC1 <- car::Anova(m.PC1, type = 3, singular.ok = T)
  eta.PC1 <- eta_squared(A.PC1, partial = F, alternative = "two.sided")
  
  m.PC2 <- lm(paste("bill.PC2 ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  A.PC2 <- car::Anova(m.PC2, type = 3, singular.ok = T)
  eta.PC2 <- eta_squared(A.PC2, partial = F, alternative = "two.sided")
  
  m.weight <- lm(paste("weight ~ ", lmvars, sep = ""), data = df.finch.sub, na.action = na.exclude)
  A.weight <- car::Anova(m.weight, type = 3, singular.ok = T)
  eta.weight <- eta_squared(A.weight, partial = F, alternative = "two.sided")
  
  # add sign to es ----------------------------------------------------------
  #get coefficients for PC1 and PC2
  df.coef <- tibble(locus.mod = names(m.PC1$coefficients),
                    coef = m.PC1$coefficients,
                    coef2 = m.PC2$coefficients,
                    coef3 = m.weight$coefficients)
  
  df.coef$locus <- str_sub(df.coef$locus.mod,1,nchar(df.coef$locus.mod)-1)
  
  df.coef.sign <- df.coef %>% 
    filter(locus!="(Intercept") %>% 
    group_by(locus) %>% 
    summarise(sum.coef = sum(coef),
              sum.coef2 = sum(coef2),
              sum.coef3 = sum(coef3)) %>% 
    mutate(sign = ifelse(sum.coef >= 0, "+",
                         ifelse(sum.coef <0, "-", NA))) %>% 
    mutate(sign2 = ifelse(sum.coef2 >= 0, "+",
                         ifelse(sum.coef2 <0, "-", NA))) %>% 
    mutate(sign3 = ifelse(sum.coef3 >= 0, "+",
                          ifelse(sum.coef3 <0, "-", NA))) %>% 
    select(-sum.coef, -sum.coef2, -sum.coef3)
  
  # three dimensions --------------------------------------------------------
  
  df.effect.size <- tibble(
    `Bill depth` = eta.depth$Eta2,
    `Bill width` = eta.width$Eta2,
    `Bill length` = eta.length$Eta2,
    locus = target_loci2
  ) %>% 
    pivot_longer(cols = -locus, values_to = "effect.size", names_to = "phenotype") %>% 
    mutate(locus = factor(locus, levels = target_loci2))
  

  p1 <- ggplot(data = df.effect.size,
               aes(label = round(effect.size, 2),x = locus, y = phenotype, fill = effect.size)) + 
    geom_tile() +
    geom_text() +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA,
                        limits = c(0,.20)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank(),
          legend.position = "none")
  
  
  
  # PCs ---------------------------------------------------------------------
  
  df.effect.sizePCs <- tibble(`Bill PC1` = eta.PC1$Eta2,
                              `Bill PC1 low` = eta.PC1$CI_low,
                              `Bill PC1 high` = eta.PC1$CI_high,
                              `Bill PC2` = eta.PC2$Eta2,
                              `Bill PC2 low` = eta.PC2$CI_low,
                              `Bill PC2 high` = eta.PC2$CI_high,
                              `Weight` = eta.weight$Eta2,
                              `Weight low` = eta.weight$CI_low,
                              `Weight high` = eta.weight$CI_high,
                              `Bill PC1 p` = A.PC1$`Pr(>F)`[2:(length(target_loci2)+1)],
                              `Bill PC2 p` = A.PC2$`Pr(>F)`[2:(length(target_loci2)+1)],
                              `Weight p` = A.weight$`Pr(>F)`[2:(length(target_loci2)+1)],
                              locus = target_loci2,
                              `Bill PC1 var` = rep(var(df.finch.sub$bill.PC1), length(target_loci2)),
                              `Bill PC2 var` = rep(var(df.finch.sub$bill.PC2), length(target_loci2))
  ) 
  #add sign
  
  df.effect.sizePCs <- df.effect.sizePCs%>% 
    pivot_longer(cols = -locus, values_to = "parameter", names_to = "phenotype") %>% 
    mutate(locus = factor(locus, levels = target_loci))
  

# add variation? ----------------------------------------------------------
  
  df.finch.sub2 <- df.finch.sub
  df.finch.sub2 <- df.finch.sub2 %>% 
    select(all_of(target_loci))
  df.finch.sub2$index <- 1:nrow(df.finch.sub2)
  
  #df.finch.sub2$genotype <- df.finch.sub2$HMGA2.simple
  
  df.finch.sub2[df.finch.sub2=="2"]<-"AB"
  df.finch.sub2[df.finch.sub2=="1"]<-"AA"
  df.finch.sub2[df.finch.sub2=="3"]<-"BB"
  
  df.finch.sub2.long <- df.finch.sub2 %>% 
    pivot_longer(cols = -index, values_to = "genotype", names_to = "locus") %>% 
    select(-index)
  
  
  df.freq <- df.finch.sub2.long %>% group_by(locus, genotype) %>% 
    summarise(count = n(), .groups = "keep")
  
  df.freq.haps <- df.freq %>% 
    pivot_wider(values_from = count, names_from = genotype) %>% 
    replace(is.na(.), 0) %>% 
    mutate(n = AA + AB + BB+ `NA`) %>% 
    mutate(af.A = ((2*AA) + AB)/(2*(n-`NA`)),
           af.B = ((2*BB) + AB)/(2*(n-`NA`)),
           q = AA / (n - `NA`), #genotype frequency of AA
           p = BB / (n - `NA`),
           af.sum = af.A + af.B,
           MAF = min(af.A, af.B),
           bAF = af.B,
           minGF = min(q, p))
  
  df.effect.sizePCs <- left_join(df.effect.sizePCs, df.freq.haps[,c("locus","MAF","bAF", "minGF")], by = "locus")

  miss_af <- df.freq.haps %>% filter(!locus %in% target_loci2) %>% select("locus","MAF","bAF", "minGF")
  # -------------------------------------------------------------------------

  #update sign iteratively
  df.effect.sizePCs <- left_join(df.effect.sizePCs, df.coef.sign, by = "locus")
  
  #update sign PC1
  df.effect.sizePCs_A <- df.effect.sizePCs %>% 
    filter(grepl("PC1", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign == "-", parameter * -1, parameter))
  
  #update sign PC2
  df.effect.sizePCs_B <- df.effect.sizePCs %>% 
    filter(grepl("PC2", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign2 == "-", parameter * -1, parameter))
  
  #update sign weight
  df.effect.sizePCs_C <- df.effect.sizePCs %>% 
    filter(grepl("Weight", phenotype) & !grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    mutate(parameter = ifelse(sign3 == "-", parameter * -1, parameter))
  
  #extract p values
  df.effect.sizePCs_p <- df.effect.sizePCs %>% 
    filter(grepl(" p", phenotype) |  grepl("var", phenotype)) 
  
  #recombine 3 dfs
  df.effect.sizePCs <- rbind(df.effect.sizePCs_A, df.effect.sizePCs_B, df.effect.sizePCs_C, df.effect.sizePCs_p) %>% 
    select(-sign, -sign2, -sign3)
  
  p2 <- df.effect.sizePCs %>% 
    filter(!grepl("high", phenotype) & !grepl("low", phenotype)) %>% 
    filter(!grepl(" p", phenotype) & !grepl("var", phenotype)) %>% 
    ggplot(aes(label = round(parameter, 2),x = locus, y = phenotype, fill = parameter)) + 
    geom_tile() +
    geom_text() +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA,
                        limits = c(0,.20)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  #p1/p2
  #ggsave(paste0("output/GEMMA/effect_size_custom_",species,"_long.pdf"), width = 8, height = 2)
  
  p2
  ggsave(paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_long.pdf"), width = 11, height = 8.5)
  
  list.plots <- list()
  
  for (genotype in target_loci){
       
    list.plots[[genotype]] <- df.finch.sub %>% 
      filter(!is.na(get(genotype))) %>% 
      ggplot(aes_string(x = genotype, y = "bill.PC1")) +
      geom_jitter(width = 0.2) + 
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      theme_bw()
  }
  
  pZ <- wrap_plots(list.plots, nrow = 3)
  pZ
  ggsave(paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_", output.name,"_boxplot.pdf"), width = 11, height = 8.5)
  
  df.effect.sizePCs.wide <- df.effect.sizePCs %>% 
    pivot_wider(values_from = parameter, names_from = phenotype) %>% 
    mutate(across(where(is.numeric), round, 3)) 
  
  #make null data for loci that were not variable enough to be included (set to 0 for pretty plotting)
  
  miss_loci <- target_loci[!target_loci %in% target_loci2]

  df_null <- tibble(MAF = miss_af$MAF,
                    bAF = miss_af$bAF,
                    minGF = miss_af$minGF,
                              `Bill PC1` = 0,
                              `Bill PC1 low` = 0,
                              `Bill PC1 high` = 0,
                              `Bill PC2` = 0,
                              `Bill PC2 low` = 0,
                              `Bill PC2 high` = 0,
                              `Weight` = 0,
                              `Weight low` = 0,
                              `Weight high` = 0,
                              `Bill PC1 p` = NA,
                              `Bill PC2 p` = NA,
                              `Weight p` = NA,
                              locus = miss_af$locus,
                              `Bill PC1 var` = rep(var(df.finch.sub$bill.PC1), nrow(miss_af)),
                              `Bill PC2 var` = rep(var(df.finch.sub$bill.PC2), nrow(miss_af))
  ) 
  
  df.effect.sizePCs.wide <- rbind(df.effect.sizePCs.wide, df_null)
  
  df.effect.sizePCs.wide$`Bill PC1 p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Bill PC1 p`, method = "BH")
  df.effect.sizePCs.wide$`Bill PC2 p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Bill PC2 p`, method = "BH")
  df.effect.sizePCs.wide$`Weight p.adjust` <- p.adjust(df.effect.sizePCs.wide$`Weight p`, method = "BH")
  
  write_csv(df.effect.sizePCs.wide, paste0("output/GEMMA/effect_sizes/PCs_effect_size_custom_",species,"_",output.name,"_es.csv"))
  
}


plot_effect_size("fortis", gwas_loci , "gwas_loci", "cluster1")
plot_effect_size("scandens", gwas_loci , "gwas_loci", "cluster2")
plot_effect_size("magnirostris", gwas_loci , "gwas_loci", "cluster3")

plot_effect_size("fortis", ld_filt_loc , "ld_prune", "cluster1")
plot_effect_size("scandens", ld_filt_loc , "ld_prune", "cluster2")
plot_effect_size("magnirostris", ld_filt_loc , "ld_prune", "cluster3")

plot_effect_size("fortis", gwas_loci , "gwas_loci", "cluster1", TRUE)
plot_effect_size("scandens", gwas_loci , "gwas_loci", "cluster2", TRUE)

plot_effect_size("fortis", gwas_loci, "gwas_loci_residuals","cluster1", resid = TRUE)
plot_effect_size("scandens", gwas_loci, "gwas_loci_residuals", "cluster2", resid = TRUE)
plot_effect_size("magnirostris", gwas_loci, "gwas_loci_residuals", "cluster3", resid = TRUE)

# -------------------------------------------------------------------------

plot_ef <- function(es_path, output.name){
  
  df.es.wide <- read.csv(es_path)
  df.es.wide <- left_join(df.es.wide, ld_filt_loc.names, by = "locus")
  df.es.wide <- df.es.wide %>% arrange(Bill.PC1) %>% mutate(order = 1:n())
  df.es.wide$tidy_locus <- fct_reorder(df.es.wide$tidy_locus, df.es.wide$order)
  
  
  pPC1 <- ggplot(df.es.wide, aes(y = tidy_locus, x = Bill.PC1,
                                 xmin = Bill.PC1.low, xmax = Bill.PC1.high)) +
    geom_errorbarh(height=.1) +
    geom_point() +
    geom_point(data = subset(df.es.wide, Bill.PC1.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) +
    labs(y = "Locus", x = "Effect size Bill PC1")
  
  #df.es.wide <- df.es.wide %>% arrange(Bill.PC2) %>% mutate(order = 1:n())
  #df.es.wide$locus <- fct_reorder(df.es.wide$locus, df.es.wide$order)
  
  pPC2 <- ggplot(df.es.wide, aes(y = tidy_locus, x = Bill.PC2,
                                 xmin = Bill.PC2.low, xmax = Bill.PC2.high)) +
    geom_errorbarh(height=.1) +
    geom_point() +
    geom_point(data = subset(df.es.wide, Bill.PC2.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) +
    labs(y = NULL, x = "Effect size Bill PC2")
  
  pWeight <- ggplot(df.es.wide, aes(y = tidy_locus, x = Weight,
                                 xmin = Weight.low, xmax = Weight.high)) +
    geom_errorbarh(height=.1) +
    geom_point() +
    geom_point(data = subset(df.es.wide, Weight.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) +
    labs(y = "Locus", x = "Effect size Weight")
  
  pPC1 + pPC2 + pWeight
  
  ggsave(paste0("output/GEMMA/effect_sizes/effect_size_", output.name,".pdf"), width = 11, height = 8.5)
  ggsave(paste0("output/GEMMA/effect_sizes/effect_size_", output.name,".png"), width = 11, height = 8.5)
  
}

plot_ef("output/GEMMA/effect_sizes/PCs_effect_size_custom_fortis_ld_prune_es.csv", "ldprune_fortis_es_noresid")
plot_ef("output/GEMMA/effect_sizes/PCs_effect_size_custom_fortis_gwas_loci_es.csv", "gwas_fortis_es_noresid")

plot_ef("output/GEMMA/effect_sizes/PCs_effect_size_custom_scandens_gwas_loci_es.csv", "scandens_gwas_loci_es_noresid")

plot_ef("output/GEMMA/effect_sizes/PCs_effect_size_custom_magnirostris_gwas_loci_es.csv", "magnirostris_gwas_loci_es_noresid")

plot_ef("output/GEMMA/effect_sizes/PCs_effect_size_custom_fuliginosa_gwas_loci_es.csv", "fuliginosa_gwas_loci_es_noresid")

# all species -------------------------------------------------------------

plot_multispecies_es <- function(df.in, output.name){
  df.es.wide <- df.in
  df.es.wide <- left_join(df.es.wide, ld_filt_loc.names, by = "locus")
  
  df.es.wide$Bill.PC1.p.adjust2 <- p.adjust(df.es.wide$Bill.PC1.p, method = "BH")
  df.es.wide$Bill.PC2.p.adjust2 <- p.adjust(df.es.wide$Bill.PC2.p, method = "BH")
  df.es.wide$Weight.p.adjust2 <- p.adjust(df.es.wide$Weight.p, method = "BH")
  
  #order by fortis
  df.loc.order <- df.es.wide %>% filter(Species == "fortis") %>% 
    arrange(Bill.PC1) %>% 
    mutate(order = 1:n()) %>% 
    dplyr::select(tidy_locus, order) %>% 
    mutate(tidy_locus = fct_reorder(tidy_locus, order))
  #right_join(df.es.wide, by = "locus")
  
  #df.es.wide$locus <- fct_reorder(df.es.wide$locus, df.es.wide$order)
  
  df.es.wide$tidy_locus <- factor(df.es.wide$tidy_locus, levels = df.loc.order$tidy_locus)
  df.es.wide$Species <- factor(df.es.wide$Species, levels = c("magnirostris","scandens", "fortis"))
  
  df.es.wide<-df.es.wide %>% mutate(Bill.PC1.sig = ifelse(Bill.PC1.p.adjust2<0.05, "Y", "N"),
                                    Bill.PC2.sig = ifelse(Bill.PC2.p.adjust2<0.05, "Y", "N"),
                                    Weight.sig = ifelse(Weight.p.adjust2<0.05, "Y", "N"))
  
  pPC1 <- ggplot(df.es.wide, aes(y = tidy_locus, x = Bill.PC1,
                                 xmin = Bill.PC1.low, xmax = Bill.PC1.high,
                                 color = Species, fill = Species, shape = Bill.PC1.sig)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    #geom_point(data = subset(df.es.wide, Bill.PC1.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none") +
    labs(y = "Locus", x = "Effect size Bill PC1") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    scale_shape_manual(values = c(19,8))
  pPC1
  
  #df.es.wide <- df.es.wide %>% arrange(Bill.PC2) %>% mutate(order = 1:n())
  #df.es.wide$locus <- fct_reorder(df.es.wide$locus, df.es.wide$order)
  
  pPC2 <- ggplot(df.es.wide, aes(y = tidy_locus, x = Bill.PC2,
                                 xmin = Bill.PC2.low, xmax = Bill.PC2.high,
                                 color = Species, fill = Species, shape = Bill.PC2.sig)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    #geom_point(data = subset(df.es.wide, Bill.PC2.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) +
    labs(y = NULL, x = "Effect size Bill PC2") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    theme(axis.text.y = element_blank(),
          legend.position = "none") +
    scale_shape_manual(values = c(19,8))
  pPC2
  
  pWeight <- ggplot(df.es.wide, aes(y = tidy_locus, x = Weight,
                                    xmin = Weight.low, xmax = Weight.high,
                                    color = Species, fill = Species, shape = Weight.sig)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    #geom_point(data = subset(df.es.wide, Bill.PC2.p.adjust < 0.05), color = "red") +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines")) +
    labs(y = NULL, x = "Effect size Weight") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    theme(axis.text.y = element_blank(),
          legend.position = "none") +
    scale_shape_manual(values = c(19,8))
  pWeight
  
  df.es.wide$Species2 <- factor(df.es.wide$Species, levels = c("fortis", "scandens","magnirostris"))
  df.es.wide$Species3 <- factor(df.es.wide$Species, levels = c("magnirostris", "scandens","fortis"))
  
  pMAF <- ggplot(df.es.wide, aes(y = tidy_locus, x = Species2,
                                 fill = MAF, label = MAF, 2)) +
    geom_tile() +
    geom_text(color = "black", size = 2.5) +
    theme_minimal() +
    theme(text=element_text(size=18, color="black")) +
    labs(y = NULL, x = "Species") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA, limits = c(0,1))
  
  pMAF.bar <- ggplot(df.es.wide, aes(y = tidy_locus, x = MAF,
                                     fill = Species3)) +
    geom_col(width = 0.55, position = position_dodge(width=0.65)) +
    theme_minimal() +
    theme(text=element_text(size=18, color="black"),
          legend.position = "none",
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c(magn.color, scan.color, fort.color)) +
    labs(y = NULL, x = "MAF") +
    scale_x_continuous(limits = c(0,.5), breaks = c(0,.5))
  
  #do genotype freqs instead of MAF
  pGenofreq <- ggplot(df.es.wide, aes(y = tidy_locus, x = Species2,
                                      fill = minGF, label = minGF, 2)) +
    geom_tile() +
    geom_text(color = "black", size = 2.5) +
    theme_minimal() +
    theme(text=element_text(size=18, color="black")) +
    labs(y = NULL, x = "Species") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(face = "italic", angle = 90)) +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA, limits = c(0,1),
                        name = "q") #not not technically q - this is minor genotype frequency within species
  pGenofreq.bar <- ggplot(df.es.wide, aes(y = tidy_locus, x = minGF,
                                          fill = Species3)) +
    geom_col(width = 0.65, position = position_dodge(1)) +
    theme_minimal() +
    theme(text=element_text(size=18, color="black"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = c(magn.color, scan.color, fort.color)) +
    labs(y = NULL, x = "q")
  
  layout <- "AAABBBCCCD"
  
  pPC1 + pPC2 + pWeight+ pMAF.bar +
    plot_layout(design = layout)
  
  ggsave(paste0("output/GEMMA/effect_sizes/effect_size_", output.name,".pdf"), width = 11, height = 6)
  ggsave(paste0("output/GEMMA/effect_sizes/effect_size_", output.name,".png"), width = 11, height = 6)
  
}

plot_multispecies_es(df.es.wide, "three_species_gwas_loci")

# clusters instead of species ---------------------------------------------

es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_es.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_gwas_loci_es.csv")
es_sca$Species <- "scandens"
es_mag <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster3_gwas_loci_es.csv")
es_mag$Species <- "magnirostris"

df.es.wide_clusters <- rbind(es_for, es_sca, es_mag) 

plot_multispecies_es(df.es.wide_clusters, "three_species_gwas_loci_clusters")

df.es.wide_clusters %>% 
  group_by(Species) %>% 
  summarise(PC1 = sum(abs(Bill.PC1)),
            #PC1_low = sum(abs(Bill.PC1.low)),
            #PC1_high = sum(abs(Bill.PC1.high)),
            PC2 = sum(abs(Bill.PC2)))
            #PC2_low = sum(abs(Bill.PC2.low)),
            #PC2_high = sum(abs(Bill.PC2.high)))

df.es.wide_clusters %>% 
  filter(locus != "HMGA2.simple" & Species == "fortis") %>% 
  select(Species, locus, Bill.PC1, Bill.PC2) %>% 
  arrange(Bill.PC1)

df.es.wide_clusters %>% 
  filter(locus == "HMGA2.simple") %>% 
  select(Species, Bill.PC1, Bill.PC2)


# plot just cluster 1 -----------------------------------------------------

es_for <- left_join(es_for, ld_filt_loc.names, by = "locus")

#order by fortis
df.loc.order <- es_for %>% filter(Species == "fortis") %>% 
  arrange(Bill.PC1) %>% 
  mutate(order = 1:n()) %>% 
  dplyr::select(tidy_locus, order) %>% 
  mutate(tidy_locus = fct_reorder(tidy_locus, order))

es_for$tidy_locus <- factor(es_for$tidy_locus, levels = df.loc.order$tidy_locus)
es_for$Species <- factor(es_for$Species, levels = c("magnirostris","scandens", "fortis"))

es_for<-es_for %>% mutate(Bill.PC1.sig = ifelse(Bill.PC1.p.adjust<0.05, "Y", "N"),
                                  Bill.PC2.sig = ifelse(Bill.PC2.p.adjust<0.05, "Y", "N"),
                                  Weight.sig = ifelse(Weight.p.adjust<0.05, "Y", "N"))


forPC1 <- es_for %>% 
  ggplot(aes(y = tidy_locus, x = Bill.PC1,
                                 xmin = Bill.PC1.low, xmax = Bill.PC1.high,
                                 color = Species, fill = Species, shape = Bill.PC1.sig)) +
  geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
  geom_point(position = position_dodge(width=0.65), size = 2) +
  #geom_point(data = subset(df.es.wide, Bill.PC1.p.adjust < 0.05), color = "red") +
  geom_vline(xintercept=0, color="black",  alpha=.5) +
  theme_minimal()+
  theme(text=element_text(size=18, color="black"))+
  theme(panel.spacing = unit(1, "lines"),
        legend.position = "none") +
  labs(y = "Locus", x = "Effect size Bill Size (PC1)") +
  scale_color_manual(values = c(fort.color)) +
  scale_shape_manual(values = c(8))

forPC2 <- es_for %>% 
  ggplot(aes(y = tidy_locus, x = Bill.PC2,
             xmin = Bill.PC2.low, xmax = Bill.PC2.high,
             color = Species, fill = Species, shape = Bill.PC2.sig)) +
  geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
  geom_point(position = position_dodge(width=0.65), size = 2) +
  #geom_point(data = subset(df.es.wide, Bill.PC1.p.adjust < 0.05), color = "red") +
  geom_vline(xintercept=0, color="black",  alpha=.5) +
  theme_minimal()+
  theme(text=element_text(size=18, color="black"))+
  theme(panel.spacing = unit(1, "lines"),
        legend.position = "none") +
  labs(y = NULL, x = "Effect size Bill Shape (PC2)") +
  scale_color_manual(values = c(fort.color)) +
  scale_shape_manual(values = c(19,8))

forPC1 + forPC2
ggsave("output/GEMMA/effect_sizes/fortis_cluster1_effect_size.png", width = 8, 
       height = 5)

# with residuals ----------------------------------------------------------

es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_residuals_es.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_gwas_loci_residuals_es.csv")
es_sca$Species <- "scandens"
es_mag <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster3_gwas_loci_residuals_es.csv")
es_mag$Species <- "magnirostris"

df.es.wide_clusters_residuals <- rbind(es_for, es_sca, es_mag) 

plot_multispecies_es(df.es.wide_clusters_residuals, "three_species_gwas_loci_clusters_residuals")

df.es.wide_clusters_residuals %>% 
  group_by(Species) %>% 
  summarise(PC1 = sum(abs(Bill.PC1)),
            PC2 = sum(abs(Bill.PC2)))

df.es.wide_clusters_residuals %>% 
  filter(Species == "fortis") %>% 
  select(Species, locus, Bill.PC1, Bill.PC2) %>% 
  arrange(Bill.PC1)


# multi species for all loci ----------------------------------------------

es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_ld_prune_es.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_ld_prune_es.csv")
es_sca$Species <- "scandens"
es_mag <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster3_ld_prune_es.csv")
es_mag$Species <- "magnirostris"

df.es.wide <- rbind(es_for, es_sca, es_mag) 

plot_multispecies_es(df.es.wide, "three_species_ld_prune_loci")

# make a legend -----------------------------------------------------------

df.es.wide2 <- rbind(es_for, es_sca, es_mag, es_ful) 
df.es.wide2$Species2 <- factor(df.es.wide2$Species, levels = c("fortis", "fuliginosa", "scandens","magnirostris"))

p.tmp <- ggplot(df.es.wide2, aes(y = locus, x = Bill.PC2, color = Species2)) +
  geom_point() +
  scale_color_manual(values = c(fort.color, fuli.color, scan.color, magn.color),
                     name = "Species", 
                     labels = c("*G. fortis*", "*G. fuliginosa*", 
                                "*G. scandens*", "*G. magnirostris*")) +
  theme(legend.position = "bottom",
        legend.text = element_markdown()) #from ggtext for markdown formatting labels

legend <- cowplot::get_legend(p.tmp)

pdf(paste0("output/GEMMA/effect_sizes/effect_size_", output.name,"_LEGEND.pdf"), width = 5, height = 2)
grid.newpage()
grid.draw(legend)
dev.off()

df.es.wide %>% 
  group_by(Species) %>% 
  summarise(PC1_sum_es = sum(abs(Bill.PC1)),
            PC2_sum_es = sum(abs(Bill.PC2)))

# is the B allele the common mag allele -----------------------------------

es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_ld_prune_es.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_ld_prune_es.csv")
es_sca$Species <- "scandens"
es_mag <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster3_ld_prune_es.csv")
es_mag$Species <- "magnirostris"
es_ful <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_fuliginosa_ld_prune_es.csv")
es_ful$Species <- "fuliginosa"

df.es.wide3 <- rbind(es_for, es_sca, es_mag, es_ful) 
df.es.wide3 <- left_join(df.es.wide3, ld_filt_loc.names, by = "locus")

df.es.wide3 %>% 
  ggplot(aes(x = tidy_locus, y = bAF, 
             fill = factor(Species, levels = c("fuliginosa", "fortis", "scandens", "magnirostris")))) +
  geom_col(position = position_dodge(), color = "white") +
  scale_fill_manual(values = c(fuli.color, fort.color, scan.color, magn.color)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="Species")) +
  labs(x = NULL, y = "Frequency of the large fortis allele")

ggsave("output/combined_genotypes/frequency_of_ld_pruned_loci_four_species_Ballele.pdf", width =8, height = 4)




# with ancestry filter ----------------------------------------------------

es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_ancfilt_es.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_gwas_loci_ancfilt_es.csv")
es_sca$Species <- "scandens"
df.es.wide_ancfilt <- rbind(es_for, es_sca) 

plot_multispecies_es(df.es.wide_ancfilt, "three_species_gwas_loci_clusters_ancfilt")

df.comp <- inner_join(df.es.wide_clusters, df.es.wide_ancfilt, by = c("locus", "Species"))
df.comp <- left_join(df.comp, ld_filt_loc.names, by = "locus")

p1<-df.comp %>% 
  filter(Species == "fortis") %>% 
  ggplot(aes(x = Bill.PC1.x, y = Bill.PC1.y, color = tidy_locus)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(y = "With ancestry filter", x = "Without ancestry filter") +
  geom_abline(intercept = 0) +
  lims(x = c(0,0.35), y = c(0,0.35)) +
  ggtitle("fortis") +
  geom_errorbar(aes(ymin = Bill.PC1.low.y,ymax = Bill.PC1.high.y)) + 
  geom_errorbarh(aes(xmin = Bill.PC1.low.x,xmax = Bill.PC1.high.x)) +
  theme(legend.title = element_blank())

p2<-df.comp %>% 
  filter(Species == "scandens") %>% 
  ggplot(aes(x = Bill.PC1.x, y = Bill.PC1.y, color = tidy_locus)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(y = "With ancestry filter", x = "Without ancestry filter") +
  geom_abline(intercept = 0) +
  lims(x = c(0,0.35), y = c(0,0.35)) +
  ggtitle("scandens") +
  geom_errorbar(aes(ymin = Bill.PC1.low.y,ymax = Bill.PC1.high.y)) + 
  geom_errorbarh(aes(xmin = Bill.PC1.low.x,xmax = Bill.PC1.high.x)) +
  theme(legend.title = element_blank())


p1+p2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("output/GEMMA/effect_sizes/effect_size_ancestry_filter_comparison.pdf",width = 8, height = 6)


# including ancestry as var -----------------------------------------------


es_for <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_with_ancestry_es_no_residuals_with_ancestry.csv")
es_for$Species <- "fortis"
es_sca <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster2_gwas_loci_with_ancestry_es_no_residuals_with_ancestry.csv")
es_sca$Species <- "scandens"
df.es.wide_ancfilt <- rbind(es_for, es_sca) 

df.comp <- inner_join(df.es.wide_clusters, df.es.wide_ancfilt, by = c("locus", "Species"))
df.comp <- left_join(df.comp, ld_filt_loc.names, by = "locus")

p3<-df.comp %>% 
  filter(Species == "fortis") %>% 
  ggplot(aes(x = Bill.PC1.x, y = Bill.PC1.y, color = tidy_locus)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(y = "With ancestry as covariate", x = "Without ancestry filter") +
  geom_abline(intercept = 0) +
  lims(x = c(0,0.35), y = c(0,0.35)) +
  ggtitle("fortis") +
  geom_errorbar(aes(ymin = Bill.PC1.low.y,ymax = Bill.PC1.high.y)) + 
  geom_errorbarh(aes(xmin = Bill.PC1.low.x,xmax = Bill.PC1.high.x)) +
  theme(legend.title = element_blank())

p4<-df.comp %>% 
  filter(Species == "scandens") %>% 
  ggplot(aes(x = Bill.PC1.x, y = Bill.PC1.y, color = tidy_locus)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(y = "With ancestry as covariate", x = "Without ancestry filter") +
  geom_abline(intercept = 0) +
  lims(x = c(0,0.35), y = c(0,0.35)) +
  ggtitle("scandens") +
  geom_errorbar(aes(ymin = Bill.PC1.low.y,ymax = Bill.PC1.high.y)) + 
  geom_errorbarh(aes(xmin = Bill.PC1.low.x,xmax = Bill.PC1.high.x)) +
  theme(legend.title = element_blank())

(p1+p2) / (p3 + p4) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
ggsave("output/GEMMA/effect_sizes/effect_size_ancestry_filter_comparison_4panel.pdf",width = 8, height = 6)
ggsave("output/GEMMA/effect_sizes/effect_size_ancestry_filter_comparison_4panel.png",width = 8, height = 6)


# compare with and without sex --------------------------------------------

es_for1 <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_residuals_nosex_es.csv")
es_for1$Species <- "fortis"

es_for2 <- read.csv("output/GEMMA/effect_sizes/PCs_effect_size_custom_cluster1_gwas_loci_residuals_es.csv")
es_for2$Species <- "fortis"

df.comp1 <- inner_join(es_for1, es_for2, by = c("locus", "Species"))
df.comp1 <- left_join(df.comp1, ld_filt_loc.names, by = "locus")

df.comp1 %>% 
  filter(Species == "fortis") %>% 
  ggplot(aes(x = Bill.PC1.x, y = Bill.PC1.y, color = tidy_locus)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(y = "With sex as covariate", x = "Without sex as covariate") +
  geom_abline(intercept = 0) +
  #lims(x = c(0,0.35), y = c(0,0.35)) +
  ggtitle("fortis") +
  geom_errorbar(aes(ymin = Bill.PC1.low.y,ymax = Bill.PC1.high.y)) + 
  geom_errorbarh(aes(xmin = Bill.PC1.low.x,xmax = Bill.PC1.high.x)) +
  theme(legend.title = element_blank())
ggsave("output/GEMMA/effect_sizes/effect_size_residuals_with_or_without_sex.png",width = 6, height = 6)