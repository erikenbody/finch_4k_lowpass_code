library(tidyverse)
library(patchwork)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

df.dap <- df.all.genos %>% filter(Island == "Daphne")

df.lgef.af <- read.csv("output/combined_genotypes/allele_frequency_annual.csv")

# make list of 1983 samples -----------------------------------------------

df.dap %>% 
  filter(Species == "fortis" & First.year.min == 1983) %>% 
  select(meaningful.unique) %>% 
  write_tsv("data/randsnps/fortis_1983.txt", col_names = F)

df.dap %>% 
  filter(Species == "fortis") %>% 
  select(meaningful.unique) %>% 
  write_tsv("data/randsnps/fortis.txt", col_names = F)


df.dap %>% filter(Species == "scandens") %>% select(meaningful.unique) %>% 
  write_tsv("data/randsnps/scandens.txt", col_names = F)

df.dap %>% 
  filter(Species == "scandens" & First.year.min == 1983) %>% 
  select(meaningful.unique) %>% 
  write_tsv("data/randsnps/scandens_1983.txt", col_names = F)

df.dap %>% filter(Species == "magnirostris") %>% select(meaningful.unique) %>% 
  write_tsv("data/randsnps/magnirostris.txt", col_names = F)

df.dap %>% 
  filter(Species == "magnirostris" & First.year.min == 1983) %>% 
  select(meaningful.unique) %>% 
  write_tsv("data/randsnps/magnirostris_1983.txt", col_names = F)

#get starting af of each locus
df.lgef.af.f <- df.lgef.af %>% 
  filter(year == 1983) %>% 
  select(Species, tidy_locus, maf) %>% 
  mutate(min_maf = maf - 0.05, max_maf = maf + 0.05)

# -------------------------------------------------------------------------

#outside R go subset VCF

# get af ------------------------------------------------------------------

#this does it from bins of 0.1 up to 0.5
getaf <- function(species){
  df.af <- read.table(paste0("data/randsnps/",species,"_1983_af.pos"), header = F)
  
  df.af <- df.af %>% 
    mutate(maf = ifelse(V5 >= 0.5, 1-V5, V5)) %>% 
    filter(maf > 0) %>% 
    mutate(points_bin = cut(maf, breaks=c(0.05, .1, .2, .3, .4, .5)))
  df.af %>% group_by(points_bin) %>% summarise(n = n())
  
  set.seed(42)
  df.af.select <- df.af %>% 
    group_by(points_bin) %>% 
    slice_sample(n = 100) 
  
  df.af.select %>%   ungroup() %>% 
    filter(!is.na(points_bin)) %>% 
    select(V1, V2) %>% 
    arrange(V1, V2) %>% 
    write_tsv(paste0("data/randsnps/",species,"_1983_selected_snps.txt"), col_names = F)
  
  df.af.bin <- df.af.select %>% 
    mutate(locus = paste(V1, V2, sep = "_")) %>% 
    select(locus, points_bin)
  
  df.af.bin <- df.af.bin %>% filter(!is.na(points_bin))
  return(df.af.bin)
}

df.af.fort <- getaf("fortis")
df.af.fort$Species <- "fortis"
df.af.scan <- getaf("scandens")
df.af.scan$Species <- "scandens"
df.af.magn <- getaf("magnirostris")
df.af.magn$Species <- "magnirostris"

df.afbins <- rbind(df.af.fort, df.af.scan, df.af.magn)


# get af from starting AF of gwas loci -------------------------------------------------

getaf2 <- function(species){
  df.af <- read.table(paste0("data/randsnps/",species,"_1983_af.pos"), header = F)
  
  df.lgef.sp <- df.lgef.af.f %>% filter(Species == species)
  df.lgef.sp.w <- df.lgef.sp %>% 
    select(Species, tidy_locus, maf) %>% 
    pivot_wider(names_from = tidy_locus, values_from = maf)
  df.af$Species <- species
  
  df.af <- left_join(df.af, df.lgef.sp.w, by = "Species")
  df.af <- df.af %>% 
    mutate(maf = ifelse(V5 >= 0.5, 1-V5, V5)) 
    
  mega.list <- list()
  for(locus in unique(df.lgef.sp$tidy_locus)){
    print(locus)
    mega.list[[locus]] <- df.af %>% 
      filter(maf < !!sym(locus) + 0.01 & maf > !!sym(locus) - 0.01) %>% 
      mutate(tidy_locus = locus)
  }
  
  mega.df <- bind_rows(mega.list)
  
  mega.df %>% group_by(tidy_locus) %>% summarise(n = n())
  
  set.seed(42)
  df.af.select <- mega.df %>% 
    group_by(tidy_locus) %>% 
    slice_sample(n = 100) 
  
  df.af.select %>% ungroup() %>% 
    filter(!is.na(tidy_locus)) %>% 
    select(V1, V2) %>% 
    arrange(V1, V2) %>% 
    distinct(.keep_all = T) %>% #some loci will be in the same points bin
    write_tsv(paste0("data/randsnps/",species,"_1983_selected_snps_refined.txt"), col_names = F)
  
  df.af.bin <- df.af.select %>% 
    mutate(locus = paste(V1, V2, sep = "_")) %>% 
    select(Species, locus, tidy_locus)
  
  df.af.bin <- df.af.bin %>% filter(!is.na(tidy_locus))
  return(df.af.bin)
}

df.af2.fort <- getaf2("fortis")
df.af2.scan <- getaf2("scandens")
#df.af.magn <- getaf2("magnirostris") #fails because of no invariant SNPs to match to the loci that are at maf =0

df.afbins2 <- rbind(df.af2.fort, df.af2.scan)

# -------------------------------------------------------------------------

#go subset the VCF to get these sites

# load back in genos ------------------------------------------------------

getgenos <- function(species){
  df.genos <- read.table(paste0("data/randsnps/",species,"_sel_snps_refined.pos.genos"), header = T) #updated to be the set of variants in af bins based on sel loci
  df.genos <- df.genos[,1:(length(df.genos)-1)]#there is a trailing NA record
  
  df.genos.long <- df.genos %>%
    pivot_longer(cols = -c("CHROM", "POS", "ALT", "REF"), names_to = "sample", values_to = "genos") %>% 
    mutate(sample = gsub("X", "", sample))
  
  df.genos.wide <- df.genos.long %>% 
    mutate(idx = paste(CHROM, POS, sep = "_")) %>% 
    select(-REF, -ALT, -CHROM, -POS) %>% 
    pivot_wider(names_from = idx, values_from = genos)
  
  df.genos.wide <- as.data.frame(df.genos.wide)
  df.genos.wide[df.genos.wide==0]<-"AA"
  df.genos.wide[df.genos.wide==1]<-"AB"
  df.genos.wide[df.genos.wide==2]<-"BB"
  return(df.genos.wide)
}

df.geno.fort <- getgenos("fortis")

df.geno.scan <- getgenos("scandens")

#df.geno.magn <- getgenos("magnirostris")

# -------------------------------------------------------------------------
meta.df <- df.all.genos %>% 
  select(meaningful.unique, First.year.min, Last.year, Species)

calcaf <- function(df.genos.wide, species){
  
  df.genos.wide <- left_join(df.genos.wide, meta.df, by = c("sample" = "meaningful.unique"))
  df.genos.wide$meaningful.unique <- df.genos.wide$sample
  
  
  df.tmp <- df.genos.wide %>% select(starts_with("chr"))
  
  these.genos <- names(df.tmp)
  
  mega.list <- list() 
  
  for(geno.names in these.genos){
  
    df.genos.wide$geno <- df.genos.wide[,geno.names]
    df.year.spread <- df.genos.wide %>% 
      filter(!is.na(First.year.min)) %>% filter(!is.na(Last.year)) %>% 
      filter(!is.na(geno)) %>% # this is a pretty strong filter but it ensures that af is calculated against the number of individuals genotyped
      select(meaningful.unique, geno, First.year.min, Last.year, Species) %>% 
      mutate(year = map2(First.year.min, Last.year, `:`)) %>% 
      select(-First.year.min, -Last.year) %>% 
      unnest(cols = c(year))
    
    df.sum.year2 <- df.year.spread %>% group_by(year, Species) %>% 
      summarise(year.n = n(), .groups = "keep") 
    
    df.year.spread$geno <- as.factor(df.year.spread$geno)
    
    df.year.prop2 <- df.year.spread %>% 
      group_by(year, geno, Species) %>% 
      summarise(geno_sum = n(), .groups = "keep") %>% 
      left_join(df.sum.year2, by = c("year", "Species")) %>% 
      mutate(prop = geno_sum / year.n)
    
    geno.match <- tibble(geno = c("AA", "AB", "BB", "BC","AC", "CC"),
                         geno.num = c(0,1,2,1,1,2))
    
    df.year.prop3 <- left_join(df.year.prop2, geno.match, by = "geno")
    
    # -calculate minor allele freq -------------------------------------------------------------
    
    df.af <- df.year.prop3 %>% 
      filter(grepl("B", geno)) %>% 
      mutate(dip.geno = geno_sum * geno.num) %>% 
      group_by(year, year.n, Species) %>% 
      summarise(geno.sum = sum(dip.geno, na.rm = T), .groups = "keep") %>% 
      mutate(af = geno.sum / (year.n * 2))
    
    df.af.1983 <- df.af %>% 
      filter(year == 1983)
    
    if(nrow(df.af.1983)>=1){
      df.af <-df.af.1983 %>% #get early time point with large number of samples
        ungroup() %>% 
        select(Species, af) %>% 
        mutate(min_allele = ifelse(af > 0.5, "A", "B")) %>% #set minor allele as the rarer 
        select(-af) %>% 
        left_join(df.af, by = "Species") %>% #add back to AF df
        mutate(maf = ifelse(min_allele == "A", 1-af, af)) #calc minor allele freq
    }else{
      df.af <- df.af %>% 
        mutate(min_allele = "B",
               maf = 1-af, af)
    }
    
    #if(minor.allele == "A"){
    #  df.af$maf <- 1-df.af$maf
    #}
    
    df.af[,"locus"] <- geno.names
    mega.list[[geno.names]] <- df.af
    
    print(geno.names)
  }
  
  df.finch.af.rand <- bind_rows(mega.list)

  return(df.finch.af.rand)
}

df.af.fort2 <- calcaf(df.geno.fort, "fortis")
df.af.scan2 <- calcaf(df.geno.scan, "scandens")
#df.af2.magn <- calcaf(df.geno.magn, "magnirostris")

df.finch.af.rand <- rbind(df.af.fort2, df.af.scan2)

#df.af.scan2 %>% filter(!locus %in% df.af2.scan$locus)

#left_join(df.af.scan2, df.af2.scan, by = "locus")

#get af bin
df.finch.af.rand <- left_join(df.finch.af.rand, df.afbins2, by = c("locus", "Species"))

write_csv(df.finch.af.rand, "output/randsnps/random_snps_refined_bined.csv")

df.finch.af.rand.sum <- df.finch.af.rand %>% 
  group_by(Species, year, tidy_locus) %>% 
  summarise(mean_af = mean(maf), .groups = "keep")
  
# plot --------------------------------------------------------------------

df.finch.af.rand %>% 
  filter(year > 1983 & Species == "fortis") %>% 
  ggplot(aes(x = year, y = maf)) +
  geom_smooth(aes(color = tidy_locus))


p.mafrand <- df.finch.af.rand.sum %>% 
  filter(year > 1983 & Species == "fortis" | year > 1983 & Species == "scandens" | year > 1983 & Species == "magnirostris") %>% 
  ggplot(aes(x = year, y = mean_af)) +
  geom_line(aes(group = tidy_locus, color = tidy_locus)) +
  geom_point(aes(group = tidy_locus)) +
  scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  ylab("Minor allele frequency") +
  facet_grid(~Species)
p.mafrand


# delta af ----------------------------------------------------------------

df.finch.af.rand.sum.f <- df.finch.af.rand.sum %>% 
  filter(year > 1982 & Species == "fortis" | year > 1982 & Species == "scandens" | year > 1982 & Species == "magnirostris") 

df.finch.af.rand.sum.f <- df.finch.af.rand.sum.f %>% 
  group_by(tidy_locus, Species) %>% 
  mutate(lag_maf = lag(mean_af))

df.finch.af.rand.sum.f$delta_maf <- df.finch.af.rand.sum.f$mean_af - df.finch.af.rand.sum.f$lag_maf

dv.avgchange.rand <- df.finch.af.rand.sum.f %>% 
  group_by(year, Species) %>% 
  summarise(mean_delta_af = mean(abs(delta_maf), na.rm = T),
            min_delta_af = min(abs(delta_maf), na.rm = T),
            max_delta_af = max(abs(delta_maf), na.rm = T)) %>% 
  filter(year > 1984) %>% 
  mutate(Species = factor(Species, levels = c("fortis", "scandens", "magnirostris")))

p.avgchange2 <-  dv.avgchange.rand %>% 
  ggplot(aes(x = year, y = mean_delta_af, color = Species)) +
  geom_line(size = 1.5, alpha = 0.6) +
  scale_x_continuous(limits = c(1982.5, 2012.5), 
                     breaks = seq(1985, 2010, by = 5)) +
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
  labs(y = expression(Mean ~ Delta ~ MAF))

p.avgchange2 + ylim(0,.1)

dv.avgchange.rand %>% 
  write_csv("output/randsnps/rand_snps_delta_af_annual_refined.csv")
