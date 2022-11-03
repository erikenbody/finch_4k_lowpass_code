library(tidyverse)
library(patchwork)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

gwas_loci <- c("gwas_genotype_chr1_2", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18", 
               "gwas_genotype_chr9_20")

#tidy_names <- c("GF01", "GF29", "GF07", "GF03", "GF30", "GF27")
tidy_names <- c("GF01", "GF03", "GF07", "GF27", "GF29", "GF30")

df.tmp <- data.frame(locus = gwas_loci)
ld_filt_loc.names <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")

df.tmp2 <- left_join(df.tmp, ld_filt_loc.names, by = "locus")



df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

df.dap <- df.all.genos %>% filter(Island == "Daphne")

head(df.dap)

#get fortis and individuals in fortis pop
df.dap.f <- df.dap %>% 
  filter(Species == "fortis" | meaningful.unique == "01Dap18974" | meaningful.unique == "00Dap19781") %>% 
  select(meaningful.unique, First.year.min, Last.year, all_of(gwas_loci))

###

df.dap.f <- df.dap.f %>%
  mutate(drought_surv = ifelse(First.year.min <= 2004 & Last.year == 2004, "died",
                               ifelse(First.year.min <= 2004 & Last.year > 2004, "survived", NA))) %>% 
  filter(!is.na(drought_surv))

# -------------------------------------------------------------------------
ftable(df.dap.f$drought_surv)


write.csv(df.dap.f, "output/drought_analysis/table_of_drought_fortis.csv")


df.dap.f.wide <- df.dap.f %>% 
  pivot_longer(cols = -c(meaningful.unique, First.year.min, 
                         Last.year, drought_surv), names_to = "locus", values_to = "genotype")

df.total_sample <- df.dap.f.wide %>% 
  group_by(locus, genotype) %>% 
  summarise(n.total = n())

df.dap.f.wide.sum <- df.dap.f.wide %>% 
  filter(!is.na(genotype)) %>% 
  mutate(genotype = as.factor(genotype), drought_surv = as.factor(drought_surv)) %>% 
  group_by(locus,genotype, drought_surv, .drop = FALSE) %>% 
  summarise(n = n())

df.dap.f.wide.sum %>% 
  filter(locus == "gwas_genotype_chr2_17")

df.dap.f.wide.sum <- left_join(df.dap.f.wide.sum, df.total_sample, by = c("locus", "genotype")) %>% 
  mutate(prop = n / n.total)

df.dap.f.wide.sum <- left_join(df.dap.f.wide.sum, ld_filt_loc.names, by = "locus")

df.dap.f.wide.sum$drought_surv <- factor(df.dap.f.wide.sum$drought_surv, levels = c("survived", "died"))

#calc sel coef
df.survived <- df.dap.f.wide.sum %>% 
  filter(genotype == "AA" | genotype == "BB") %>% 
  filter(drought_surv == "survived") %>% 
  arrange(tidy_locus, genotype)

df.sel <- df.survived %>% 
  group_by(tidy_locus) %>% 
  summarise(fitness = prop / sum(prop),
            sel_coef = 1 - fitness) 

df.total <- df.dap.f.wide %>% 
  group_by(locus) %>% 
  summarise(n.full = n())

df.sel2 <- left_join(df.dap.f.wide.sum, df.total, by = "locus")


####USE THIS##############
df.Vll <- df.sel2 %>% 
  select(-prop) %>% 
  filter(genotype!="AB") %>% 
  pivot_wider(names_from = drought_surv, values_from = n) %>% 
  mutate(GF_pre = (died + survived) / n.full,
         GF_post = (survived)/n.full) %>% 
  select(locus, genotype, GF_pre, GF_post) %>% 
  pivot_wider(values_from = c("GF_pre", "GF_post"), 
              names_from = "genotype") %>% 
  mutate(Vll = (GF_post_BB * GF_pre_AA) / (GF_post_AA * GF_pre_BB),
         sel = 1- Vll)

df.Vll <- left_join(df.Vll, ld_filt_loc.names, by = "locus")
df.Vll %>%select(sel, tidy_locus) %>%  arrange(tidy_locus)

###########################

# fisher test -------------------------------------------------------------

df.fet <- data.frame(pval = rep(NA, 6))
row.names(df.fet) <- unique(df.sel2$tidy_locus)
for(focus in unique(df.sel2$tidy_locus)){
  
  df.tmp <- df.sel2 %>% filter(tidy_locus == focus)
  print(focus)
  df.tmp <- df.tmp %>% 
    ungroup() %>% 
    select(drought_surv, genotype, n) %>% 
    mutate(A = ifelse(genotype == "AA", 2*n,
                      ifelse(genotype == "BB", 0,
                             ifelse(genotype == "AB", 1*n, NA))),
           B = ifelse(genotype == "AA", 0,
                       ifelse(genotype == "BB", 2*n,
                              ifelse(genotype == "AB", 1*n, NA)))) %>% 
    group_by(drought_surv) %>% 
    summarise(Atot = sum(A), Btot = sum(B)) %>% 
    data.frame()
  
  row.names(df.tmp) <- df.tmp$drought_surv ; df.tmp$drought_surv <- NULL
  print(df.tmp)
  fet <- fisher.test(df.tmp)
  print(fet)
  df.fet[focus, "pval"] <- round(fet$p.value, 4)
  df.fet[focus, "estimate"] <- round(fet$estimate, 4)
  
}



# -------------------------------------------------------------------------




df.survived$sel_coef <- df.sel$sel_coef
df.AA.sel <- df.survived %>% 
  filter(genotype == "BB")


df.dap.f.wide.sum %>%   
  ggplot(aes(x = drought_surv, y = prop, color = genotype, group = genotype)) + 
  geom_point(position = position_dodge(0.04)) +
  geom_line(position = position_dodge(0.04)) +
  geom_text(data = df.AA.sel, aes(x = 1.5, y = 0.8, label = round(sel_coef, 2)),
            color = "black") +
  facet_wrap(~tidy_locus) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        text = element_text(size = 16)) + 
  labs(x = NULL, y = "Proportion of individuals")
ggsave("output/drought_analysis/survival_geno_plot.png", width = 8.5, height = 4)
ggsave("output/drought_analysis/survival_geno_plot.pdf", width = 8.5, height = 4)



# dot matrix --------------------------------------------------------------


df.dap.f.wide <- left_join(df.dap.f.wide, ld_filt_loc.names, by = "locus")

for(loc in tidy_names){
  print(loc)
  df.tmp <- df.dap.f.wide %>% filter(tidy_locus == !!(loc))
  
  plots <- list()
  
  for(geno in c("AA","AB","BB")){
    df.tmp.geno <- df.tmp %>% filter(genotype == !!geno)
    start_lat <- 0
    start_lng <- 0
    griddf <- expand.grid(x1 = seq(from = start_lat, by = 1, l = nrow(df.tmp.geno)/4),
                          y1 = seq(from = start_lng, by = 1, l = nrow(df.tmp.geno)/4))
    
    df.tmp.geno <- df.tmp.geno %>% arrange(drought_surv)
    df.tmp.geno$x1 <- griddf[1:nrow(df.tmp.geno),]$x1
    df.tmp.geno$y1 <- griddf[1:nrow(df.tmp.geno),]$y1
    
    plots[[geno]] <- df.tmp.geno %>% 
      ggplot() + 
      geom_point(aes(x = factor(x1), y = factor(y1), color = drought_surv), size = 25) +
      theme_void() +
      theme(legend.position = "none") +
      scale_color_manual(values = c("#f04f4f", "grey70")) #+
      #ggtitle(paste0(loc,"_",geno)) 
    #ggsave(paste0("output/drought_analysis/dotmats/",loc,"_",geno,".pdf"), width = 6, height = 6)
  }

  if(loc == "GF30"){
    design <- "AAAAA
             BBBB#
             CC###"
    wrap_plots(plots, nrow = 3, design = design)
    ggsave(paste0("output/drought_analysis/dotmats/",loc,".pdf"), width = 12, height = 12)
  } else{
    design <- "AA###
             BBBBB
             CC###"
    wrap_plots(plots, nrow = 3, design = design)
    ggsave(paste0("output/drought_analysis/dotmats/",loc,".pdf"), width = 12, height = 12)
  }
  

}

start_lat <- 0
start_lng <- 0
griddf <- expand.grid(x1 = seq(from = start_lat, by = 1, l = 10),
                      y1 = seq(from = start_lng, by = 1, l = 8))
plot(griddf)

df.tmp <- df.dap.f.wide %>% filter(locus == "HMGA2.simple") %>% 
  arrange(genotype)
df.tmp$x1 <- griddf[1:nrow(df.tmp),]$x1
df.tmp$y1 <- griddf[1:nrow(df.tmp),]$y1

df.tmp %>% 
  ggplot() + 
  geom_point(aes(x = x1, y = y1, color = genotype), size = 3) +
  theme_void()


df.dap.f.wide %>% 
  #filter(locus == "HMGA2.simple") %>% 
  select(meaningful.unique, drought_surv, genotype, locus) %>% 
  filter(!is.na(genotype)) %>% 
  group_by(genotype) %>% 
  arrange(drought_surv) %>% 
  mutate(ind = 1:n()) %>% 
  ggplot() +
  geom_point(aes(x = ind, y = genotype, color = drought_surv), size = 2,
              width = 0.2, height=0.2) +
  theme_void() +
  theme(axis.text.x = element_blank()) +
  facet_grid(~locus)
ggsave("output/drought_analysis/dotmatrix_in.pdf", width = 11, height = 8.5)

# barplot -----------------------------------------------------------------

df.dap.f.wide.sum %>%   
  filter(drought_surv == "survived") %>% 
  ggplot(aes(x = genotype, y = prop, fill = genotype)) + 
  geom_col(position = position_dodge(), color = "black") + 
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none") + 
  facet_wrap(~tidy_locus) +
  labs(y = "Proportion survived", x = NULL) +
  scale_fill_manual(values = c("#fee8c8", "#fdbb84","#e34a33"))
ggsave("output/drought_analysis/barchart_six_loci.pdf", width = 7.5, height = 4)

library(ggpattern)
df.dap.f.wide.sum %>%   
  filter(drought_surv == "survived") %>% 
  mutate(genotype = factor(gsub("B","L",gsub("A","S", genotype)), levels = c("SS","SL","LL"))) %>% 
  ggplot(aes(x = genotype, y = prop, fill = genotype)) + 
  geom_col_pattern(
    aes(pattern = genotype, pattern_angle = genotype, pattern_spacing = genotype), 
    fill            = 'white',
    colour          = 'black', 
    pattern_density = 0.1, 
    pattern_fill    = 'black',
    pattern_colour  = 'black'
  ) +
  theme_bw() +
  scale_pattern_spacing_discrete(range = c(0.04, 0.05)) + 
  theme(legend.position = 'none') + 
  #coord_fixed(ratio = 1) + 
  facet_wrap(~tidy_locus) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none",
        text = element_text(size = 18)) + 
  labs(y = "Proportion survived", x = NULL) 

ggsave("output/drought_analysis/barchart_six_loci_patterns.pdf", width = 7.5, height = 4)


# or just plain -----------------------------------------------------------

df.dap.f.wide.sum %>%   
  filter(drought_surv == "survived") %>% 
  mutate(genotype = factor(gsub("B","L",gsub("A","S", genotype)), levels = c("SS","SL","LL"))) %>% 
  ggplot(aes(x = genotype, y = prop, fill = genotype)) + 
  geom_col(fill = "grey80", color = "black") +
  theme_bw() +
  theme(legend.position = 'none') + 
  facet_wrap(~tidy_locus) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none",
        text = element_text(size = 18)) + 
  labs(y = "Proportion survived", x = NULL) 

ggsave("output/drought_analysis/barchart_six_loci_grey.pdf", width = 7.5, height = 4)

tmp<-df.dap.f.wide.sum %>% filter(drought_surv == "survived" & locus == "HMGA2.simple" )
ggplot(tmp, aes(x="", y=prop, fill=genotype)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void()



# sum of all allleles -----------------------------------------------------
df.dap.f.num<-df.dap.f
df.dap.f.num[df.dap.f.num=="AB"]<-1
df.dap.f.num[df.dap.f.num=="AA"]<-2
df.dap.f.num[df.dap.f.num=="BB"]<-0

df.num.loc <- df.dap.f.num %>% select(all_of(gwas_loci))
df.num.loc <- mutate_all(df.num.loc, function(x) as.numeric(as.character(x)))
df.dap.f.num$sum.alt <- rowSums(df.num.loc, na.rm = T)

df.dap.f.num <- df.dap.f.num %>% 
  mutate(survival = as.numeric(gsub("survived", 1, gsub("died", 0, drought_surv)))) 
df.dap.f.num %>% 
  ggplot(aes(x = sum.alt, y = survival)) + 
  geom_jitter(height = 0.01, width = 0) +
  stat_smooth(method="glm", se=T, method.args = list(family=binomial)) +
  labs(x = "Genotype dosage", y = "Survival") +
  theme_bw()
ggsave("output/drought_analysis/survival_curve_genotype_dosage.png", width = 8, height = 6)

t.test(df.dap.f.num$sum.alt ~ df.dap.f.num$drought_surv)

full_mod_AIC <- AIC(glm(survival ~ sum.alt, data = df.dap.f.num, family = "binomial"))

#try dropping loci
for(loc in names(df.num.loc)){
  df.tmp <- df.num.loc[,-which(names(df.num.loc) %in% c(loc))]
  df.dap.f.num[,paste0("drop_",loc)] <- rowSums(df.tmp, na.rm = T)
  partial_AIC <- AIC(glm(survival ~ unlist(df.dap.f.num[,paste0("drop_",loc)]), data = df.dap.f.num, family = "binomial"))
  if(partial_AIC < full_mod_AIC){
    print(paste0("drop_",loc))
    print(partial_AIC)
  }
}

#update without G07 and G27
df.dap.f.num$sum.alt.filt <- rowSums(df.num.loc[,-which(names(df.num.loc) %in% 
                                                          c("ALX1.simple", "gwas_genotype_chr9_20"))], na.rm = T)


part_mod_AIC <- AIC(glm(survival ~ sum.alt.filt, data = df.dap.f.num, family = "binomial"))


g03_only_AIC <- AIC(glm(survival ~ HMGA2.simple, data = df.dap.f.num, family = "binomial"))

# boxplot -----------------------------------------------------------------
#filtered
df.dap.f.num %>% 
  ggplot(aes(x = drought_surv, y = sum.alt.filt)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.4,
               binwidth = .05, color = "grey") +
  #geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = NULL, y = "Number of small alleles") +
  theme_bw() +
  theme(text = element_text(size = 18)) 
  
ggsave("output/drought_analysis/survival_boxplot_genotype_dosage.png", width = 3, height = 6)
ggsave("output/drought_analysis/survival_boxplot_genotype_dosage.pdf", width = 3, height = 6)


# -------------------------------------------------------------------------
#all 6 loci

df.dap.f.num %>% 
  ggplot(aes(x = drought_surv, y = sum.alt)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.4,
               binwidth = .2, color = "grey") +
  #geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = NULL, y = "Number of small alleles") +
  theme_bw() +
  theme(text = element_text(size = 18)) 
