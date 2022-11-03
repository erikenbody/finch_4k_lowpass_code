library(tidyverse)

path1 <- "data/gcta/greml_PC1_maf05.txt.hsq"
qtlpath <- "data/gcta/greml_PC1_with_qtls_maf05.txt.hsq"

calc_vqtl <- function(path1, qtlpath){
  
  df1 <- read_table(path1)
  df1$type <- "base"
  dfQ <- read_table(qtlpath)
  dfQ$type <- "QTL"
  
  df.comp <- rbind(df1, dfQ)
  df.comp <- df.comp %>% 
    filter(!is.na(SE)) %>% 
    pivot_wider(names_from = type, values_from = c(Variance, SE))
  
  df.comp %>% 
    filter(Source == "V(G)") %>% 
    mutate(Vqtl = Variance_base - Variance_QTL,
           diffSE = SE_base - SE_QTL,
           total_var = Vqtl + Variance_QTL + mean(SE_base, SE_QTL)) %>% 
    return()
  
}

#warnings are OK

calc_vqtl("data/gcta/greml_PC1_maf05.txt.hsq", "data/gcta/greml_PC1_with_qtls_maf05.txt.hsq")
calc_vqtl("data/gcta/greml_PC2_maf05.txt.hsq", "data/gcta/greml_PC2_with_qtls_maf05.txt.hsq")
calc_vqtl("data/gcta/greml_weight_INts_maf05.txt.hsq", "data/gcta/greml_weight_INts_with_qtls_maf05.txt.hsq")
