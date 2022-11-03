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


gwas_loci <- c("gwas_genotype_chr1_1", "gwas_genotype_chr1A_19", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_21", 
               "gwas_genotype_chr9_23")



df.sub <- df.dap.genos %>% 
  filter(Species == "fortis") %>% 
  select(First.year.min, Last.year,contains("anc"), all_of(gwas_loci))

df.sub <-df.sub %>% 
  mutate(across(gwas_loci, ~ifelse( .x == "AB", 1,
                                    ifelse(.x == "BB", 2,
                                           ifelse(.x == "AA", 0, NA)))))

df.sub %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(HMGA2.simple), y = anc_fuli)) +
  theme_bw() 




df.sub <- df.sub %>% 
  mutate(bin_geno = ifelse(HMGA2.simple == "AB" | HMGA2.simple == "BB", 1, 
                           ifelse(HMGA2.simple == "AA", 0, NA))) %>% 
  filter(!is.na(bin_geno))

df.sub %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(bin_geno), y = anc_fuli)) +
  theme_bw() +
  labs(x = "Carries at least 1 large HMGA2 allele", y = "G. fuliginosa ancestry")

t.test(anc_fuli ~ bin_geno, data = df.sub)


# scan --------------------------------------------------------------------

df.sub <- df.dap.genos %>%
  filter(Species == "scandens") %>% 
  select(contains("anc"), ALX1.simple)

df.sub <- df.sub %>% 
  mutate(bin_geno = ifelse(ALX1.simple == "AB" | ALX1.simple == "BB", 1, 
                           ifelse(ALX1.simple == "AA", 0, NA))) %>% 
  filter(!is.na(bin_geno))

df.sub %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(bin_geno), y = anc_fort)) +
  theme_bw() +
  labs(x = "Carries at least 1 blunt ALX1 allele", y = "G. fortis ancestry")

t.test(anc_fort ~ bin_geno, data = df.sub)

df.tmp <- df.all.genos %>% 
  filter(grm_cluster == "cluster1")
cor.test(df.tmp$bill.PC1.sp, df.tmp$weight)
