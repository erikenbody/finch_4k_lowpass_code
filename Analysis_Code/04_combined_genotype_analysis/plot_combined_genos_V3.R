library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggsci)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")
head(df.all.genos)
df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")


# -------------------------------------------------------------------------
#gwas_loci <- c("gwas_genotype_chr1_1", "gwas_genotype_chr1A_19", "ALX1.simple",
#               "HMGA2.simple", "gwas_genotype_chr2_21", 
#               "gwas_genotype_chr9_23")

#V2
gwas_loci <- c("gwas_genotype_chr1_2_ext", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18_ext", 
               "gwas_genotype_chr9_20")

#hand curated by LD profile and highest pvalue SNP in fortis gwas
ld_filt_loc <- c(gwas_loci, "asm28loci_chr2_10", "asm28loci_chr2_12","asm28loci_chr2_17",
                 "asm28loci_chr2_19","asm28loci_chr3_22", "asm28loci_chr3_24", "asm28loci_chr5_25",
                 "asm28loci_chr7_26", "asm28loci_chr25_28")


df.tmp <- df.dap.genos %>% dplyr::select(all_of(ld_filt_loc))

these.genos <- names(df.tmp)

calc_af <- function(these.genos, maf_year, minor_allele){
  
  #maf year toggles between minor allele as starting af in 198X, or 
  #minor allele relative to small allele
  
  mega.list <- list() 
  
  for(geno.names in these.genos){
    
    df.dap.genos$geno <- df.dap.genos[,geno.names]
    df.year.spread <- df.dap.genos %>% 
      filter(!is.na(First.year.min)) %>% filter(!is.na(Last.year)) %>% 
      filter(!is.na(geno)) %>% # this is a pretty strong filter but it ensures that af is calculated against the number of individuals genotyped
      dplyr::select(meaningful.unique, geno, First.year.min, Last.year, Species) %>% 
      mutate(year = map2(First.year.min, Last.year, `:`)) %>% 
      dplyr::select(-First.year.min, -Last.year) %>% 
      unnest(cols = c(year))
    
    df.sum.year2 <- df.year.spread %>% group_by(year, Species) %>% 
      summarise(year.n = n(), .groups = "keep") 
    
    df.year.prop2 <- df.year.spread %>% 
      group_by(year, geno, Species) %>% 
      summarise(geno_sum = n(), .groups = "keep") %>% 
      left_join(df.sum.year2, by = c("year", "Species")) %>% 
      mutate(prop = geno_sum / year.n)
    
    geno.match <- tibble(geno = c("AA", "AB", "BB", "BC","AC", "CC"),
                         geno.num = c(0,1,2,1,1,2))
    
    df.year.prop3 <- left_join(df.year.prop2, geno.match, by = "geno")
    
    
    # identify minor allele early in study ------------------------------------
    #some franken code to force MAF to fortis
    #BB.freq <- df.year.prop3 %>% filter(Species == "fortis" & year == 1988 & geno == "BB")
    #BA.freq <- df.year.prop3 %>% filter(Species == "fortis" & year == 1988 & geno == "AB")
    #B.sum <- sum(max(BB.freq[1,"geno_sum"], na.rm = F) * 2, max(BA.freq[1,"geno_sum"], na.rm = T), na.rm = T)
    #B.n <- max(BA.freq[1,"year.n"], BB.freq[1,"year.n"], na.rm = T)
    #
    ##set B if it is not present at all i.e. minor allele is NA
    #minor.allele <- ifelse(B.sum / (B.n*2) > 0.5, "A", "B")
    #if(is.na(minor.allele)){
    #  minor.allele <- "B"
    #}
    
    # -calculate minor allele freq -------------------------------------------------------------
    
    df.af <- df.year.prop3 %>% 
      filter(grepl("B", geno)) %>% #when maf_year = F, this means the minor allele is the "small allele"
      mutate(dip.geno = geno_sum * geno.num) %>% 
      group_by(year, year.n, Species) %>% 
      summarise(geno.sum = sum(dip.geno, na.rm = T), .groups = "keep") %>% 
      mutate(af = geno.sum / (year.n * 2))
    
    if(maf_year == TRUE){
      df.af <- df.af %>% 
        filter(year == 1983) %>% #get early time point with large number of samples
        ungroup() %>% 
        select(Species, af) %>% 
        mutate(min_allele = ifelse(af > 0.5, "A", "B")) %>% #set minor allele as the rarer 
        select(-af) %>% 
        left_join(df.af, by = "Species") %>% #add back to AF df
        mutate(maf = ifelse(min_allele == "A", 1-af, af)) #calc minor allele freq
    } else if(maf_year == FALSE) {
      df.af <- df.af %>% 
        mutate(maf = af)
      
      if(minor_allele == "A"){
        df.af$maf <- 1-df.af$maf
      }
    }
    
    df.af[,"locus"] <- geno.names
    mega.list[[geno.names]] <- df.af
    
    print(geno.names)
  }
  
  df.finch.af <- bind_rows(mega.list)
  
  return(df.finch.af)
}

#false runs it in relation to hard coded minor allele
df.finch.af <- calc_af(these.genos, FALSE, "A")


ld_filt_loc.names <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")
df.finch.af <- left_join(df.finch.af, ld_filt_loc.names, by = "locus")

unique(df.finch.af$locus)

df.finch.af <- df.finch.af %>% 
  mutate(type = ifelse(grepl("gwas_", locus), "gwas_fortis", 
                       ifelse(grepl("asm28", locus), "asm_28_loci",
                              ifelse(grepl("chr5", locus), "chr5",
                                     ifelse(grepl("chr29", locus), "chr29", 
                                            ifelse(locus == "HMGA2.simple", "HMGA2",
                                                   ifelse(locus == "ALX1.simple", "ALX1", NA)))))))

# -------------------------------------------------------------------------

p.maf4 <- df.finch.af %>% 
  #filter(type == "asm_28_loci") %>% 
  filter(year > 1987 & Species == "fortis" | year > 1987 & Species == "scandens") %>% 
  ggplot(aes(x = year, y = maf)) +
  geom_line(aes(group = locus, color = type), size = 1, alpha = 0.6) +
  scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  ylab("Large allele frequency") + #will need to change label if switching to minor af
  facet_grid(~Species)
p.maf4
ggsave("output/combined_genotypes/annual_plots/asm28_annual_af_plots.png", width = 12, height = 8)


# stepwise ----------------------------------------------------------------

df.finch.af.f <- df.finch.af %>% 
  filter(locus %in% gwas_loci) %>% 
  filter(year > 1982 & Species == "fortis" | year > 1982 & Species == "scandens" | year > 1982 & Species == "magnirostris") 

df.finch.af.f <- df.finch.af.f %>% 
  #filter(locus == "gwas_genotype_chr1A_83") %>% 
  group_by(locus, Species) %>% 
  mutate(lag_maf = lag(maf))

df.finch.af.f$delta_maf <- df.finch.af.f$maf - df.finch.af.f$lag_maf


# get starting AF ---------------------------------------------------------

af_start <- df.finch.af %>% filter(year == 1983) %>% 
  filter(locus %in% gwas_loci) %>% 
  filter(Species == "fortis" | Species == "scandens" | Species == "magnirostris") %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens", "magnirostris"))) 


# make patchwork layout ---------------------------------------------------

layout <- "
AAAAAA
BBBBBB
BBBBBB
BBBBBB
"
# write out af ------------------------------------------------------------

write_csv(df.finch.af.f, "output/combined_genotypes/allele_frequency_annual.csv")
write_csv(df.finch.af, "output/combined_genotypes/annual_plots/all_taxa_annual_af_plots.csv")

# average changes ---------------------------------------------------------

#get random snps
df.rand <- read.csv("output/randsnps/rand_snps_delta_af_annual.csv")
df.rand <- df.rand %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens", "magnirostris")))

dv.avgchange <- df.finch.af.f %>% 
  group_by(year, Species) %>% 
  summarise(mean_delta_af = mean(abs(delta_maf), na.rm = T),
            min_delta_af = min(abs(delta_maf), na.rm = T),
            max_delta_af = max(abs(delta_maf), na.rm = T)) %>% 
  filter(year > 1984) %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens", "magnirostris")))


# calc magnitude of AF change ---------------------------------------------

df_years_change <- df.finch.af.f %>% 
  filter(year == 1983 | year == 2012) %>% 
  group_by(Species, tidy_locus, year) %>% 
  summarise(maf = maf) %>% 
  pivot_wider(values_from = maf, names_from = year) %>% 
  mutate(delta_af = round(`2012` - `1983`, 2) * 100) %>% 
  filter(Species!="magnirostris")

df_years_change %>% 
  ggplot() + geom_col(aes(x = tidy_locus, y = abs(delta_af), fill = Species), 
                      position = "dodge", width = 0.5) +
  theme_bw() 

# only fortis and scandens ------------------------------------------------
df.rand.2sp <- df.rand %>% 
  filter(Species == "fortis" | Species == "scandens") 

p.avgchange.2sp <-  dv.avgchange %>% 
  filter(Species == "fortis" | Species == "scandens") %>% 
  ggplot(aes(x = year, y = mean_delta_af, color = Species)) +
  geom_line(size = 1.5, alpha = 0.6) +
  geom_line(data = df.rand.2sp, color = "grey", size = 1.5) +
  scale_x_continuous(limits = c(1982.5, 2014.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Annual change in AF") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color)) +
  labs(y = expression( ~ Delta ~ MAF / year)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

af_start.2sp <- af_start %>% filter(Species == "fortis" | Species == "scandens")

p.maf_L2.2sp <- df.finch.af %>% 
  filter(locus %in% gwas_loci) %>% 
  filter(year > 1982 & Species == "fortis" | year > 1982 & Species == "scandens") %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens"))) %>% 
  ggplot(aes(x = year, y = maf)) +
  geom_line(aes(group = tidy_locus, color = tidy_locus), size = 1.5, alpha = 0.8) +
  #geom_point(aes(group = tidy_locus, color = tidy_locus), size = 2.5) +
  scale_x_continuous(limits = c(1982.5, 2014.5), breaks = seq(1985, 2010, by = 5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  geom_hline(data = af_start.2sp, aes(yintercept = maf, group = tidy_locus, 
                                      color = tidy_locus), linetype='dotted') +
  #geom_text(data = df_years_change, aes(x = 2014, y = `2012`, 
  #                                      label = paste0(delta_af, "%"), color = tidy_locus)) +
  geom_text_repel(data = df_years_change, aes(x = 2014, y = `2012`, 
                                        label = paste0(delta_af, "%"), 
                                        color = tidy_locus), direction = "y",
                  segment.color = "transparent", show.legend = FALSE,
                  fontface = "bold") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  #ylab("Minor allele frequency") +
  ylab("Small allele frequency") + #will need to change label if switching to minor af
  facet_grid(~Species) +
  scale_color_aaas() +
  #scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494')) +
  labs(color = "Locus") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 
#p.maf_L2.2sp

p.avgchange.2sp / p.maf_L2.2sp + plot_layout(design = layout)
ggsave("output/combined_genotypes/annual_plots/two_af_annual_plots_2sp.pdf", 
       width = 8.5, height = 7)


# -------------------------------------------------------------------------


dv.avgchange %>% 
  filter(Species == "fortis") %>% 
  ggplot(aes(x = year, y = mean_delta_af, color = Species)) +
  geom_line(size = 1.5, alpha = 0.6) +
  geom_line(data = subset(df.rand.2sp, Species == "fortis"), color = "grey", size = 1.5) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Annual change in AF") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color)) +
  labs(y = expression( ~ Delta ~ MAF / year)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave("output/combined_genotypes/just_delta_af.pdf", width = 5, height = 3)

# locus by locus --------------------------------------------------------------------

p.daf.hmga2 <- df.finch.af.f %>% 
  filter(year > 1984) %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens", "magnirostris"))) %>% 
  filter(locus == "HMGA2.simple") %>% 
  ggplot(aes(x = year, y = delta_maf)) +
  geom_line(aes(group = locus, color = Species), size = 1, alpha = 0.6) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Annual change in allele frequency") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color))

p.daf.hmga2

ggsave("output/combined_genotypes/annual_plots/hmga2_annual_change_in_af_bw.png", 
       width = 8.5, height = 6)
