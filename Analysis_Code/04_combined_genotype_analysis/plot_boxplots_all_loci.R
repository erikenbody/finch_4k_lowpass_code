library(tidyverse)

#setup data
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

#update sex with field sex if genetic is missing
df.all.genos <- df.all.genos %>% 
  mutate(genetic.sex = as.factor(gsub("unassigned", NA, genetic.sex)),
         Sex2 = ifelse(Sex == 1, "Male", 
                       ifelse(Sex ==2, "Female", NA))) %>% 
  mutate(munch.sex = ifelse(is.na(genetic.sex) & !is.na(Sex2), as.character(Sex2), as.character(genetic.sex)))


df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")

#gwas_loci <- c("gwas_genotype_chr1_1", "gwas_genotype_chr1A_19", "ALX1.simple",
#               "HMGA2.simple", "gwas_genotype_chr2_21", "gwas_genotype_chr9_23")

#V2
gwas_loci <- c("gwas_genotype_chr1_2", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18", 
               "gwas_genotype_chr9_20")


locus <- "HMGA2.simple"
species <- "fortis"

df.dap.genos[df.dap.genos=="AB"]<-"SL"
df.dap.genos[df.dap.genos=="AA"]<-"SS"
df.dap.genos[df.dap.genos=="BB"]<-"LL"

for(species in c("fortis", "scandens", "magnirostris", "fuliginosa")){
  for(locus in gwas_loci){
    for(pheno in c("bill.PC1.sp","bill.PC2.sp","weight")){
      df.dap.genos[,locus] <- factor(df.dap.genos[,locus], levels = c("SS", "SL", "LL"))
      df.dap.genos %>% 
        filter(Species == species) %>% 
        filter(!is.na(get(locus))) %>% 
        ggplot(aes(x = !!sym(locus), y = !!sym(pheno))) + 
        #geom_dotplot(binaxis = "y", stackdir = "center", 
        #             alpha = 0.5, color = "grey", dotsize = 0.3) +
        #geom_point()
        geom_jitter(color = "grey", width=0.15, height = 0) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        labs(x = NULL, y = pheno, title = species) +
        theme_bw() +
        theme(text = element_text(size = 20)) 
      ggsave(paste0("output/combined_genotypes/boxplots/",species,"_",locus,"_",pheno,".pdf"), 
             width = 5.5, height = 6)
    }
  }
}

