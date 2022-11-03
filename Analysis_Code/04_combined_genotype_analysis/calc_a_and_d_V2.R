library(tidyverse)

#setup data. 
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls_tidy.csv", na.strings = "")

#drop outliers (two fortis immigrants, one may be maginrostris)
df.all.genos <- df.all.genos %>% 
  filter(meaningful.unique!="01Dap21269" & meaningful.unique!="01Dap21274")


#update sex with field sex if genetic is missing
df.all.genos <- df.all.genos %>% 
  mutate(genetic.sex = as.factor(gsub("unassigned", NA, genetic.sex)),
         Sex2 = ifelse(Sex == 1, "Male", 
                       ifelse(Sex ==2, "Female", NA))) %>% 
  mutate(munch.sex = ifelse(is.na(genetic.sex) & !is.na(Sex2), as.character(Sex2), as.character(genetic.sex)))


# pca ---------------------------------------------------------------------

df.clust <- df.all.genos %>% 
  filter(Island == "Daphne") %>% 
  filter(!is.na(grm_cluster))

#pivot longer to get column of measurements
df.clust.long <- df.clust %>% 
  select(meaningful.unique, grm_cluster, bill.depth, bill.width, bill.length, bill.PC1.sp, bill.PC2.sp, starts_with("G")) %>% 
  select(-Genus, -genetic.sex, -Genotype) %>% 
  pivot_longer(cols = -c(meaningful.unique, starts_with("G"), grm_cluster), values_to = "value", names_to = "pheno")

#a second pivot longer to add in genotypes
df.clust.long <- df.clust.long %>% 
  pivot_longer(cols = -c(meaningful.unique, grm_cluster, pheno, value), values_to = "geno", names_to = "locus")

df.clust.calc <- df.clust.long %>% 
  filter(locus %in% c("G01", "G03", "G07", "G29", "G30", "G27")) %>% 
  group_by(pheno, grm_cluster, locus, geno) %>% 
  summarise(m = mean(value, na.rm = T),
            sd = sd(value, na.rm = T),
            n = n()) %>% 
  filter(!is.na(geno)) %>% 
  pivot_wider(names_from = geno, values_from = c("m", "sd","n")) %>% 
  mutate(a = (m_AA - m_BB)/2,
         se_a = sqrt((sd_AA/n_AA) + (sd_BB/n_BB)), #se is the square root of the sum of SD, taking into account the sample size of each group
         d = (m_AB - ((m_AA + m_BB)/2)),
         se_d = sqrt((sd_AB/n_AB) + ((sd_AA + sd_BB)/2) / ((n_AA + n_BB)/2) )) %>% #this is rough, but calculated the mean sd of AA and BB combined, as well as average sample size
  mutate(abs_diff = abs(a) - abs(d)) %>% 
  select(a, se_a, d, se_d) %>% 
  arrange(grm_cluster, pheno) 

write_csv(df.clust.calc, "output/combined_genotypes/a_d_statistics_3clusters.csv")


# plot --------------------------------------------------------------------

plot_multispecies_ad <- function(df.in, output.name){
  df.es.wide <- df.in %>% 
    filter(pheno == "bill.PC1.sp" | pheno == "bill.PC2.sp") %>% 
    mutate(tidy_locus = locus)
  
  matchdf <- tibble(grm_cluster = c("cluster1", "cluster2", "cluster3"),
         Species = c("fortis", "scandens", "magnirostris"))
  
  df.es.wide <- left_join(df.es.wide, matchdf, by = "grm_cluster")
  
  #order by fortis
  df.loc.order <- df.es.wide %>% 
    filter(grm_cluster == "cluster1" & pheno == "bill.PC1.sp") %>% 
    arrange(a) %>% 
    mutate(order = 1:n()) %>% 
    dplyr::select(tidy_locus, order) %>% 
    mutate(tidy_locus = fct_reorder(tidy_locus, order))
  #right_join(df.es.wide, by = "locus")
  
  #nix that, order by locus name:
  df.loc.order$tidy_locus <- factor(df.loc.order$tidy_locus, levels = 
                                      c("G01", "G03", "G07","G27", "G29", "G30"))
  
  #df.es.wide$locus <- fct_reorder(df.es.wide$locus, df.es.wide$order)
  
  df.es.wide$tidy_locus <- factor(df.es.wide$tidy_locus, levels =  rev(c("G01", "G03", "G07","G27", "G29", "G30")))
  df.es.wide$Species <- factor(df.es.wide$Species, levels = c("magnirostris","scandens", "fortis"))
  
  
  pPC1.a <- df.es.wide %>% filter(pheno == "bill.PC1.sp") %>% 
      ggplot(aes(y = tidy_locus, x = a,
                                   xmin = a - se_a, xmax = a + se_a,
                                   color = Species, fill = Species)) +
        geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
        geom_point(position = position_dodge(width=0.65), size = 2) +
        geom_vline(xintercept=0, color="black",  alpha=.5) +
        theme_minimal()+
        theme(text=element_text(size=18, color="black"))+
        theme(panel.spacing = unit(1, "lines"),
              legend.position = "bottom",
              axis.text.y = element_text(face = "italic"),
              axis.title.x = element_text(face = "italic"),
              legend.text = element_text(face = "italic")) +
        labs(y = "Locus", x = "a") +
        scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    lims(x = c(-2, 2)) +
    ggtitle("Beak size (PC1)")
  
  
  pPC1.d <- df.es.wide %>% filter(pheno == "bill.PC1.sp") %>% 
    ggplot(aes(y = tidy_locus, x = d,
               xmin = d - se_d, xmax = d + se_d,
               color = Species, fill = Species)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.x = element_text(face = "italic"),
          legend.text = element_text(face = "italic")) +
    labs(y = NULL, x = "d") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    lims(x = c(-2, 2))
  
  pPC1.a + pPC1.d
  ggsave(paste0("output/combined_genotypes/a_d_PC1_", output.name,".png"), width = 11, height = 6)
  
  pPC2.a <- df.es.wide %>% filter(pheno == "bill.PC2.sp") %>% 
    ggplot(aes(y = tidy_locus, x = a,
               xmin = a - se_a, xmax = a + se_a,
               color = Species, fill = Species)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "bottom",
          axis.text.y = element_text(face = "italic"),
          axis.title.x = element_text(face = "italic"),
          legend.text = element_text(face = "italic")) +
    labs(y = "Locus", x = "a") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    lims(x = c(-2, 2)) +
    ggtitle("Beak size (PC2)")
   
  
  pPC2.d <- df.es.wide %>% filter(pheno == "bill.PC2.sp") %>% 
    ggplot(aes(y = tidy_locus, x = d,
               xmin = d - se_d, xmax = d + se_d,
               color = Species, fill = Species)) +
    geom_errorbarh(height=.1, position = position_dodge(width=0.65)) +
    geom_point(position = position_dodge(width=0.65), size = 2) +
    geom_vline(xintercept=0, color="black",  alpha=.5) +
    theme_minimal()+
    theme(text=element_text(size=18, color="black"),
          axis.title.x = element_text(face = "italic"),
          legend.text = element_text(face = "italic"))+
    theme(panel.spacing = unit(1, "lines"),
          legend.position = "none",
          axis.text.y = element_blank()) +
    labs(y = NULL, x = "d") +
    scale_color_manual(values = c(magn.color, scan.color, fort.color)) +
    lims(x = c(-2, 2))
  
  pPC2.a + pPC2.d
  
  ggsave(paste0("output/combined_genotypes/a_d_PC2_", output.name,".png"), width = 11, height = 6)
  
  (  pPC1.a + pPC1.d ) /
    (pPC2.a + pPC2.d) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  ggsave(paste0("output/combined_genotypes/a_d_PC1_PC2_", output.name,".png"), width = 8, height = 10)
  
  
}
plot_multispecies_ad(df.clust.calc, "plots")
