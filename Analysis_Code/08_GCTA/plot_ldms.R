library(tidyverse)

path = "data/gcta/ldms_reml_results_2mafbins_2median_ldbins_V2.txt.hsq"
name = "None"

load_ld_hsq <- function(path, name){
  
  df.ldms <- read.table(path, fill = T, skip = 1)
  
  df.ldms1 <- df.ldms %>% filter(grepl("/Vp", V1))
  df.ldms2 <- df.ldms[12,c(1,2)]
  names(df.ldms2) <- c("V2", "V3")
  df.ldms.c <- bind_rows(df.ldms1, df.ldms2)
  
  df.ldms.c <- df.ldms.c %>% 
    mutate(V1 = ifelse(is.na(V1), "Total", V1),
           Variance = as.numeric(V2),
           SE = as.numeric(V3)) %>% 
    select(-V2,-V3)
  
  df.ldms.c$LD <- factor(c("Low LD", "High LD", "Low LD", "High LD", "Total"), levels = c("Low LD", "High LD", "Total"))
  df.ldms.c$MAF <- factor(c("Common (>0.05)", "Common (>0.05)", "Rare (<0.05)", "Rare (<0.05)", "Sum V(G)/Vp"), levels = c("Rare (<0.05)","Common (>0.05)", "Sum V(G)/Vp"))
  df.ldms.c$category <- name
  return(df.ldms.c)
}

df.ldms.c <- load_ld_hsq("data/gcta/ldms_reml_results_2mafbins_2median_ldbins_V2.txt.hsq", "None")
df.ldmsPCA <- load_ld_hsq("data/gcta/ldms_reml_results_2mafbins_2median_ldbins_20PCcov_V2.txt.hsq", "20 PCAs")
df.ldmsPCA.merge <- rbind(df.ldms.c, df.ldmsPCA)

ggplot(data = df.ldmsPCA.merge, aes(x = MAF, y = Variance, fill = category)) +
  geom_col(position = position_dodge(.9), width = 0.75) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=.2,
                position=position_dodge(.9)) +
  facet_grid(~LD, drop = TRUE, scales = "free_x") +
  theme_bw() +
  labs(y = "Variance ± SE")
ggsave("output/gcta/variance_explained_maf_ld_bins_PCA_covariates_with_grmalg1.png", width = 8, height =6 )


# plot multiple phenotypes ------------------------------------------------
df.ldms.1 <- load_ld_hsq("data/gcta/ldms_reml_results_2mafbins_2median_ldbins_V2.txt.hsq", "Bill size")
df.ldms.2 <- load_ld_hsq("data/gcta/bill_PC2_ldms_reml_results_2mafbins_2median_ldbins_V2.txt.hsq", "Bill shape")
df.ldms.3 <- load_ld_hsq("data/gcta/weightINT_ldms_reml_results_2mafbins_2median_ldbins_V2.txt.hsq", "Weight")

df.3ldms <- rbind(df.ldms.1, df.ldms.2, df.ldms.3)
df.3ldms$category <- factor(df.3ldms$category, levels = c("Weight", "Bill shape", "Bill size"))
ggplot(data = df.3ldms, aes(x = MAF, y = Variance, fill = category)) +
  geom_col(position = position_dodge(.9), width = 0.75) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=.2,
                position=position_dodge(.9)) +
  facet_grid(~LD, drop = TRUE, scales = "free_x") +
  theme_bw() +
  labs(y = "Variance ± SE") +
  theme(legend.position = "bottom")
ggsave("output/gcta/Daphne_cluster1_variance_explained_maf_ld_bins_PCA_covariates_with_grmalg1_3phenos.png", width = 8, height =6 )

# check where low LD variants come from -----------------------------------

source("code/finch_code/lowpass/15_GEMMA/gemma_functions_2021.R")
library(data.table)
gcta_lr_loco <- gcta.order("data/gcta/geno_assoc_loco.loco.mlma", "p")
gcta_lr_loco$index <- paste(gcta_lr_loco$chr, gcta_lr_loco$ps, sep = "_")

lds_seg = read.table("data/gcta/greml_ldms.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
lds_seg$chr <- gsub("31", "1A",
                     gsub("32", "4A",
                          gsub("33", "Z", lds_seg$chr)))
lds_seg$chr <- paste0("chr",lds_seg$chr)
lds_seg$index <- paste(lds_seg$chr, lds_seg$bp, sep = "_")

quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[3])
lb2 = which(lds_seg$ldscore_SNP > quartiles[3])
lb1_snp = lds_seg$index[lb1]
lb2_snp = lds_seg$index[lb2]


gcta_lr_loco <- gcta_lr_loco %>% 
  mutate(ld = ifelse(index %in% lb1_snp, "low",
                     ifelse(index %in% lb2_snp, "high", NA)))

table(gcta_lr_loco$ld)

gcta_lr_loco %>% 
  filter(!is.na(ld)) %>% 
  ggplot() +
  geom_boxplot(aes(x = ld, y = log_p))
gcta_lr_loco %>% 
  filter(chr == "chr1A") %>% 
  ggplot() + geom_point(aes(x = ps, y = log_p, color = factor(ld))) 

gcta_lr_loco %>% 
  filter(chr == "chr1") %>% 
  ggplot() + geom_point(aes(x = ps, y = log_p, color = factor(ld))) 


gcta_lr_loco %>% 
  filter(!is.na(ld)) %>% 
  ggplot() + 
  geom_point(aes(x = row, y = log_p, color = factor(ld))) +
  theme_bw() +
  theme(axis.text.x = element_blank())
ggsave("output/gcta/billPC1_manhattan_with_ldgroups_labeled.png", width = 16, height = 5)


lds_seg %>% 
  mutate(index = 1:n()) %>% 
  filter(chr == "chr4") %>% 
  ggplot(aes(x = index, y = chr, fill = ldscore_region)) + 
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "red")

