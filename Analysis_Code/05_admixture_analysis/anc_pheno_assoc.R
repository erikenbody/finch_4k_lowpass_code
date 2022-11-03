library(tidyverse)
library(broom)
#setup data
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

#get tidy names
df_tidy <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")

df.all.genos %>% 
  filter(grepl("big", Species)) %>% 
  select(ALX1.geno) %>% 
  filter(grepl("P", ALX1.geno))

# species ----------------------------------------------------------------

#df.fort4lm <- df.all.genos %>% filter(Species == "fortis" & Island == "Daphne")
df.scan4lm <- df.all.genos %>% filter(grm_cluster == "cluster1" & Island == "Daphne")

fort_pc1 <- tidy(lm(bill.PC1.sp ~ anc_scan + anc_fuli, data = df.fort4lm))
fort_pc1$pheno <- "PC1"
fort_pc1$species <- "fortis"
summary(lm(bill.PC1.sp ~ anc_scan + anc_fuli, data = df.fort4lm))$adj.r.squared

fort_pc2 <- tidy(lm(bill.PC2.sp ~ anc_scan + anc_fuli, data = df.fort4lm))
fort_pc2$pheno <- "PC2"
fort_pc2$species <- "fortis"
summary(lm(bill.PC2.sp ~ anc_scan + anc_fuli, data = df.fort4lm))$adj.r.squared

#df.scan4lm <- df.all.genos %>% filter(Species == "scandens" & Island == "Daphne")
df.scan4lm <- df.all.genos %>% filter(grm_cluster == "cluster2" & Island == "Daphne")

scan_pc1 <- tidy(lm(bill.PC1.sp ~ anc_fort + anc_fuli, data = df.scan4lm))
summary(lm(bill.PC1.sp ~ anc_fort + anc_fuli, data = df.scan4lm))
scan_pc1$pheno <- "PC1"
scan_pc1$species <- "scandens"
summary(lm(bill.PC1.sp ~ anc_fort + anc_fuli, data = df.scan4lm))$adj.r.squared

scan_pc2 <- tidy(lm(bill.PC2.sp ~ anc_fort + anc_fuli, data = df.scan4lm))
scan_pc2$pheno <- "PC2"
scan_pc2$species <- "scandens"
summary(lm(bill.PC2.sp ~ anc_fort + anc_fuli, data = df.scan4lm))$adj.r.squared

df_alltests <- rbind(fort_pc1, fort_pc2, scan_pc1, scan_pc2)
write_csv(df_alltests, "output/admixture/ancestry_phenotype_associations.csv")


# plots -------------------------------------------------------------------
library(ggpubr)

df.fort4lm.w <- df.fort4lm %>% 
  select(bill.PC1.sp, bill.PC2.sp, meaningful.unique, anc_scan, anc_fuli) %>% 
  pivot_longer(cols = -c(meaningful.unique, bill.PC1.sp, bill.PC2.sp),
               names_to = "Species", values_to = "ancestry") %>% 
  mutate(`Ancestry source` = gsub("anc_fuli","fuliginosa", 
                                  gsub("anc_scan", "scandens", Species)))

anc1 <- df.fort4lm.w %>% filter(Species == "anc_scan")
anc2 <- df.fort4lm.w %>% filter(Species == "anc_fuli")

mix_plot <- function(titlel = "G. fortis", ylabel = "Bill Size (PC1)", 
                     color1 = scan.color, color2 = fuli.color,
                     sp1 = "G. fuliginosa", sp2 = "G. scandens",
                     ypos1 = 1.3, ypos2 = -1.3, pheno = "bill.PC1.sp"){
  mix_plot <- ggplot(data = anc1, aes(x = ancestry, y = !!sym(pheno)),) + 
    geom_point(data = anc1, aes(x = ancestry, y = !!sym(pheno), color = `Ancestry source`),
               alpha = 0.8) + 
    geom_smooth(data = anc1, aes(x = ancestry, y = !!sym(pheno), color = `Ancestry source`),
                method = "lm", se=FALSE) +
    geom_point(data = anc2, aes(x = ancestry, y = !!sym(pheno), color = `Ancestry source`),
               alpha = 0.8) + 
    geom_smooth(data = anc2, aes(x = ancestry, y = !!sym(pheno), color = `Ancestry source`),
                method = "lm", se=FALSE) +
    stat_regline_equation(data = anc1, label.y = ypos1, label.x = .65, aes(label = ..eq.label..,
                                                                           color = `Ancestry source`), show.legend = FALSE) +
    stat_regline_equation(data = anc1, label.y = ypos1 - 0.3, label.x = .65, aes(label = ..rr.label..,
                                                                                 color = `Ancestry source`), show.legend = FALSE) +
    stat_regline_equation(data = anc2, label.y = ypos2 + 0.3, label.x = .65, aes(label = ..eq.label..,
                                                                                 color = `Ancestry source`), show.legend = FALSE) +
    stat_regline_equation(data = anc2, label.y = ypos2, label.x = .65, aes(label = ..rr.label..,
                                                                           color = `Ancestry source`), show.legend = FALSE) +
    theme_bw() +
    labs(x = "Ancestry", y = ylabel, title = titlel) +
    scale_color_manual(values = c(color1, color2)) +
    theme(legend.position = "bottom") +
    xlim(0,0.85)
  return(mix_plot)
}

p.fort1 <- mix_plot(color1 = fuli.color, color2 = scan.color, pheno = "bill.PC1.sp")

p.fort2 <- mix_plot(color1 = fuli.color, color2 = scan.color, pheno = "bill.PC2.sp", 
         ylabel = "Bill Shape (PC2)")

p.fort1

df.scan4lm.w <- df.scan4lm %>% 
  select(bill.PC1.sp, bill.PC2.sp, meaningful.unique, anc_fort, anc_fuli) %>% 
  pivot_longer(cols = -c(meaningful.unique, bill.PC1.sp, bill.PC2.sp),
               names_to = "Species", values_to = "ancestry") %>% 
  mutate(`Ancestry source` = gsub("anc_fuli","fuliginosa", 
                                  gsub("anc_fort", "fortis", Species)))

anc1 <- df.scan4lm.w %>% filter(Species == "anc_fort")
anc2 <- df.scan4lm.w %>% filter(Species == "anc_fuli")

p.scan1 <- mix_plot(color2 = fuli.color, color1 = fort.color, pheno = "bill.PC1.sp",
         titlel = "G. scandens", ypos1 = 1.5)

p.scan2 <- mix_plot(color2 = fuli.color, color1 = fort.color, pheno = "bill.PC2.sp",
         titlel = "G. scandens", ylabel = "Bill Shape (PC2)")

(p.fort1 + p.scan1) / (p.fort2 + p.scan2) + plot_annotation(tag_levels = 'A')
ggsave("output/admixture/ancestry_phenotype_plots.png", width = 12, height = 10)

# -------------------------------------------------------------------------

df.full2 <- read.csv("output/GEMMA/phenotype_plots/Daphne_phenotype_dataframe_for_plotting.csv")

df.dap <- df.all.genos %>% filter(Island == "Daphne") %>% select(-PC1, -PC2)
df.dap <- left_join(df.dap, df.full2[,c("meaningful.unique", "PC1", "PC2")], by = "meaningful.unique")


summary(lm(PC1 ~ anc_scan + anc_fuli + anc_fort + anc_mag + anc_mag2, data = df.dap))
