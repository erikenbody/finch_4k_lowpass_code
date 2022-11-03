library(tidyverse)
library("PerformanceAnalytics")

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

#get admix
fort_admix <- read.csv("output/admixture/fortis_admixture_longform_annual.csv")
scan_admix <- read.csv("output/admixture/scandens_admixture_longform_annual.csv")

df.admix <- rbind(fort_admix, scan_admix)

df.admix.wide <-  df.admix %>% 
  select(-n.ancestry) %>% 
  pivot_wider(values_from = "mean.ancestry", names_from = "ancestral.pop")

#get af
df.finch.af <- read.csv("output/combined_genotypes/annual_plots/all_taxa_annual_af_plots.csv")

df.finch.af.wide <- df.finch.af %>% 
  select(year, Species, maf, tidy_locus) %>% 
  pivot_wider(names_from = tidy_locus, values_from = maf)

#get phenotypes
df.pheno <- read.csv("output/phenotypic_data/annual_plots/annual_phenotype_summaries.csv")

df.merge1 <- left_join(df.admix.wide, df.finch.af.wide, by = c("year" = "year", "nest.Breeding.pop.ext" = "Species"))
df.merge2 <- left_join(df.merge1, df.pheno, by = c("year", "nest.Breeding.pop.ext" = "Species"))

df.merge2 <- df.merge2 %>% filter(year > 1982)

# correlations ------------------------------------------------------------
df.data <- df.merge2 %>% select(G01, G03, G07, G27, G29, G30,
                                fortis_mean, magnirostris_mean, scandens_mean, fuliginosa_mean,
                                bill.PC1.sp_mean, bill.PC2.sp_mean, weight_mean, nest.Breeding.pop.ext)
df.fort.sum <- df.data %>% filter(nest.Breeding.pop.ext == "fortis") %>% select(-nest.Breeding.pop.ext)
chart.Correlation(df.fort.sum, histogram=TRUE, pch=19)

cor(df.fort.sum) %>% data.frame() %>% write.csv("output/admixture/fortis_annual_correlations_ancestry_af.csv")

df.scan.sum <- df.data %>% filter(nest.Breeding.pop.ext == "scandens") %>% select(-nest.Breeding.pop.ext)
cor(df.scan.sum) %>% data.frame() %>% write.csv("output/admixture/scandens_annual_correlations_ancestry_af.csv")


# -------------------------------------------------------------------------




df.merge2 %>% 
  filter(nest.Breeding.pop.ext == "scandens") %>% 
  ggplot(aes(x = fortis_mean, y = G03)) + 
  geom_point() +
  geom_smooth(method = lm, se = F) +
  theme_bw()
  
df.merge2 %>% 
  filter(nest.Breeding.pop.ext == "fortis") %>% 
  ggplot(aes(x = fuliginosa_mean, y = G03)) + 
  geom_point() +
  geom_smooth(method = lm, se = F) +
  theme_bw() 

df.merge2 %>% 
  filter(nest.Breeding.pop.ext == "fortis") %>% 
  ggplot(aes(x = G03, y = fuliginosa_mean)) + 
  geom_point() +
  geom_smooth(method = lm, se = F) +
  theme_bw()


# predicting phenotype ----------------------------------------------------
df.merge2.fort <- df.merge2 %>% filter(nest.Breeding.pop.ext == "fortis")

cor.test(df.merge2.fort$G03, df.merge2.fort$fuliginosa_mean)
cor.test(df.merge2.fort$G07, df.merge2.fort$bill.PC1.sp_mean)
cor.test(df.merge2.fort$G03, df.merge2.fort$bill.PC1.sp_mean)
cor.test(df.merge2.fort$bill.PC2.sp_mean, df.merge2.fort$scandens_mean)

lm1 <- lm(bill.PC1.sp_mean ~ scandens_mean + fuliginosa_mean + G01 + G03 + G07 + G27 + G29 + G30, data = df.merge2.fort)
library(car)
Anova(lm1)

AIC(lm1)
AIC(lm(bill.PC1.sp_mean ~ scandens_mean + fuliginosa_mean, data = df.merge2.fort))


lm4 <- lm(bill.PC2.sp_mean ~ scandens_mean + fuliginosa_mean + G01 + G03 + G07 + G27 + G29 + G30, data = df.merge2.fort)
summary(lm4)

plot(df.merge2.fort$bill.PC1.sp_mean ~ df.merge2.fort$G07)

tidy(lm4)


# -------------------------------------------------------------------------

df.merge2.scan <- df.merge2 %>% filter(nest.Breeding.pop.ext == "scandens")
lm2 <- lm(bill.PC1.sp_mean ~ fortis_mean + fuliginosa_mean + G01 + G03 + G07 + G27 + G29 + G30, data = df.merge2.scan)
summary(lm2)

lm3 <- lm(bill.PC2.sp_mean ~ fortis_mean + fuliginosa_mean + G01 + G03 + G07 + G27 + G29 + G30, data = df.merge2.scan)
summary(lm3)

plot(df.merge2.scan$bill.PC1.sp_mean ~ df.merge2.scan$fuliginosa_mean)
plot(df.merge2.scan$bill.PC2.sp_mean ~ df.merge2.scan$fortis_mean)


