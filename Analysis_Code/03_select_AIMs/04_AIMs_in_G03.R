library(tidyverse)
library(patchwork)


fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

#setup data
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

ld_filt_loc.names <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")

#update sex with field sex if genetic is missing
df.all.genos <- df.all.genos %>% 
  mutate(genetic.sex = as.factor(gsub("unassigned", NA, genetic.sex)),
         Sex2 = ifelse(Sex == 1, "Male", 
                       ifelse(Sex ==2, "Female", NA))) %>% 
  mutate(munch.sex = ifelse(is.na(genetic.sex) & !is.na(Sex2), as.character(Sex2), as.character(genetic.sex)))


df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")


# get early ---------------------------------------------------------------
#generate list of samples to caculate Fst 

#filter(Species == species & First.year.min <= year1 & Last.year > year2)

early_samps <- df.dap.genos %>% 
  #select(locus, tidy_locus, Species, year, af) %>% 
  filter(First.year.min <= 1983 & Last.year >=1983 & Species == "scandens" |
           First.year.min <= 1983 & Last.year >=1983 & Species == "fortis" | 
           First.year.min <= 1992 & Last.year >=1992 & Species == "fuliginosa" |
           First.year.min <= 1992 & Last.year >=1992 &  Species == "magnirostris")

gwas_loci <- c("gwas_genotype_chr1_2_ext", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18_ext", 
               "gwas_genotype_chr9_20")

early_samps %>% 
  select(meaningful.unique) %>% 
  write_tsv("data/pixy/haplotype_groups/early_samples.txt", col_names = F)

for(focal in gwas_loci){
  df <- early_samps %>% 
    select(meaningful.unique, Species, !!sym(focal)) %>% 
    filter(!!sym(focal) == "AA") %>% 
    select(meaningful.unique, Species)
  
  print(focal)
  dir.create("data/pixy/haplotype_groups", showWarnings = F)
  write_tsv(df, paste0("data/pixy/haplotype_groups/", focal, "_4pixy.txt"), col_names = F)
}


# run in pixy separately --------------------------------------------------

# -------------------------------------------------------------------------
#identify SNPs with Fst > 0.5  (G03) or > 0.3 (G01) in 1983/1991

df_g3 <- read_tsv("data/pixy/four_species_pixy_G03_persite_fst.txt")
df_g1 <- read_tsv("data/pixy/four_species_pixy_G01_persite_fst.txt")
hist(df_g1$avg_wc_fst)

df_g3 %>% 
  filter(pop1 == "fortis" & pop2 == "fuliginosa") %>% 
  ggplot() +
  geom_histogram(aes(x = avg_wc_fst))

df_g1 %>% 
  filter(pop1 == "fortis" & pop2 == "fuliginosa") %>% 
  ggplot() +
  geom_histogram(aes(x = avg_wc_fst))

head(df_g3)
d1<-df_g3 %>% 
  filter(pop1 == "fortis" & pop2 == "fuliginosa") %>% 
  filter(avg_wc_fst > 0.5) %>% 
  select(chromosome, window_pos_1)
d2<-df_g1 %>% 
  filter(pop1 == "fortis" & pop2 == "fuliginosa") %>% 
  filter(avg_wc_fst > 0.3) %>% 
  select(chromosome, window_pos_1) 

#not really enough snps in G01 to be a useful analysis. Will focus on G03 from here on
nrow(d2)

rbind(d1,d2) %>% 
  write_tsv("data/pixy/G03_G01_early_pos.txt", col_names = F)

df_g3 %>% 
  filter(pop1 == "fortis" & pop2 == "fuliginosa") %>% 
  filter(avg_wc_fst > 0.5) %>% 
  summarise(minfst = min(avg_wc_fst),
            maxfst = max(avg_wc_fst))

# -------------------------------------------------------------------------

unested_years <- df.dap.genos %>%
  dplyr::filter(!is.na(First.year.min), !is.na(Last.year)) %>%
  dplyr::select(meaningful.unique, First.year.min, Last.year, Species) %>%
  dplyr::mutate(year = map2(First.year.min, Last.year, `:`)) %>%
  select(-First.year.min, -Last.year) %>%
  unnest(cols = c(year))

#set up sample list
for (fyear in c(1983:2012)){
  
  for (species in c("fortis", "scandens")){
    
    print(species)
    
    fort_unested_years <- unested_years %>% 
      dplyr::filter(Species == species)
    dir.create(paste0("data/pixy/", species, "_yearly"))
    
    fort_unested_years %>% 
      filter(year == fyear) %>% 
      select(meaningful.unique) %>% 
      write_tsv(paste0("data/pixy/", species, "_yearly/",fyear, "_", species, "_samps.txt"), col_names = F)
  }
  
}

# -------------------------------------------------------------------------
af_df <- read.table("data/pixy/fort_all_years_af.txt")
names(af_df) <- c("chromosome", "position", "ref", "alt", "allele_freq", "year")

af_df$position <- as.factor(af_df$position)

early_maf <- af_df %>% 
  filter(year == 1983 & allele_freq > 0.5)

af_df <- af_df %>% 
  mutate(allele_freq = ifelse(position %in% early_maf$position, 1-allele_freq, allele_freq))

af_df %>% 
  filter(chromosome == "chr1A" & year > 1983) %>% 
  #mutate(allele_freq = ifelse(allele_freq >=0.5, 1-allele_freq, allele_freq)) %>%
  ggplot(aes(x = year, y = allele_freq, color = position)) +
  geom_line() +
  labs(x = "Year", y = "Allele Frequency", color = "SNP Position") +
  theme_minimal()
ggsave("output/pixy/fortis_early_daf_snps_G03.png", width = 8, height = 6)

af_daf <- af_df %>% 
  filter(year == 1983 | year == 2012) %>% 
  select(chromosome, position, allele_freq, year) %>% 
  pivot_wider(names_from = year, values_from = allele_freq) %>% 
  mutate(daf = `2012` - `1983`) %>% 
  mutate(dir = ifelse(daf > 0, "+", "-")) %>% 
  group_by(chromosome) %>% 
  arrange(chromosome, -daf)

df_long <- af_daf %>%
  pivot_longer(cols = c(`1983`, `2012`), names_to = "year", values_to = "value")

# Plot the data
ggplot(df_long, aes(x = year, y = value, group = position, color = dir)) +
  geom_line() +
  geom_point() +
  #scale_color_continuous(low = "blue", high = "red") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "maf", color = "direction of change") +
  facet_grid(~chromosome)
ggsave("output/pixy/fortis_early_daf_snps_G03_binary.png", width = 8, height = 6)


# scandens ----------------------------------------------------------------

g03_annual <- af_df %>% 
  filter(chromosome == "chr1A" & year > 1982) %>% 
  #mutate(allele_freq = ifelse(allele_freq >=0.5, 1-allele_freq, allele_freq)) %>%
  ggplot(aes(x = year, y = allele_freq, group = position)) +
  geom_line() +
  labs(x = NULL, y = "Minor Allele Frequency", color = "SNP Position") +
  theme_minimal() +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5))

af_daf <- af_df %>% 
  filter(year == 1983 | year == 2012) %>% 
  select(chromosome, position, allele_freq, year) %>% 
  pivot_wider(names_from = year, values_from = allele_freq) %>% 
  mutate(daf = `2012` - `1983`) %>% 
  mutate(dir = ifelse(daf > 0, "+", "-")) %>% 
  group_by(chromosome) %>% 
  arrange(chromosome, -daf)

df_long <- af_daf %>%
  pivot_longer(cols = c(`1983`, `2012`), names_to = "year", values_to = "value")

#just plot G03
start_end_af <- df_long %>%
  filter(chromosome == "chr1A") %>% 
  ggplot(aes(x = year, y = value, group = position, color = dir)) +
  geom_line() +
  geom_point() +
  #scale_color_continuous(low = "blue", high = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Minor Allele Frequency", color = "direction of change") 

layout <- "
AAAB
"
g03_annual + start_end_af  +  plot_layout(design = layout)
ggsave("output/pixy/fortis_early_daf_snps_G03.png", width = 8, height = 6)

