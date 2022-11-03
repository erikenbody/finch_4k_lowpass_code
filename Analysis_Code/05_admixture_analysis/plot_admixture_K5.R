library(tidyverse)
library(patchwork)
library(readxl)
library(dendextend)
library(ggtext)

#df.ash.admix <- read.table("data/admixture/Daphne_FortisFuliginosaScandens_K3_AdmixProportions.txt", header = T)
df.ash.admix_og <- read.table("data/admixture/Ash_data/k5/projAdmix.5.Q", header = F)

names(df.ash.admix_og) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
sample.names <- read.table("data/admixture/Ash_data/k5/sample_order.txt", header = F)
df.ash.admix_og$sample.ID <- sample.names$V1

df.ash.admix <- read.table("data/admixture/Ash_data/admix_proportions_PRUNED/K5_proportions.txt", header = T)
df.ash.admix$sample.ID <- df.ash.admix$meaningful.unique ; df.ash.admix$meaningful.unique <- NULL

#code to get metadata from the databse
V2.mdb <- read_excel("/Users/erikenbody/Dropbox/Shared_Uppsala/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
df.final.proj <- V2.mdb %>% filter(!is.na(Lowcov_project) | !is.na(High_coverage_project))
df.final.proj %>% filter(Genus == "Geospiza" & meaningful.unique!="02Gen1389") %>%
  select(meaningful.unique, Species,Island, Genotype, nest.Breeding.pop, contains("bill"), First.year.min, Last.year) %>%
  distinct(meaningful.unique, .keep_all = T) %>%
  write_csv("data/snmf/geospiza_population_info.csv")

pop_info <- df.final.proj %>% filter(Genus == "Geospiza" & meaningful.unique!="02Gen1389") %>%
  #select(meaningful.unique) %>%
  select(meaningful.unique, Species,Island, Genotype, nest.Breeding.pop, contains("bill"), First.year.min, Last.year) %>%
  distinct(meaningful.unique, .keep_all = T)

pop_info <- read.csv("data/snmf/geospiza_population_info.csv", header = T)

df.ash.admix.long <- df.ash.admix %>% 
  pivot_longer(cols = -sample.ID)

df.ash.admix.long <- left_join(df.ash.admix.long, pop_info, by = c("sample.ID" = "meaningful.unique"))

df.ash.admix.long %>% 
  group_by(Species, name) %>% 
  summarise(mean.val = mean(value)) %>% 
  print(n = 100)


# compare og to new -------------------------------------------------------

df_comp <- inner_join(df.ash.admix_og, df.ash.admix, by ="sample.ID")

p1<-ggplot(df_comp, aes(x = Cluster3.x, y = Cluster1.y)) + geom_point() +
  ggtitle("fortis") +
  labs(x = "pruned", y= "full")

p2<-ggplot(df_comp, aes(x = Cluster2.x, y = Cluster3.y)) + geom_point() +
  ggtitle("scandens") +
  labs(x = "pruned", y= "full")

p3<-ggplot(df_comp, aes(x = Cluster4.x, y = Cluster4.y)) + geom_point() +
  ggtitle("fuliginosa") +
  labs(x = "pruned", y= "full")

p4<-ggplot(df_comp, aes(x = Cluster5.x, y = Cluster2.y)) + geom_point() +
  ggtitle("magnirostris") +
  labs(x = "pruned", y= "full")

(p1 + p2) / (p3 + p4)

# pheno -------------------------------------------------------------------

df.pheno.adx <- df.ash.admix.long %>% 
  filter(name == "Cluster3") %>% 
  filter(Species == "fortis" | Species == "scandens" | Species == "hybrid") %>% 
  select(Species, bill.length, bill.depth, value)

ggplot(data = df.pheno.adx) + geom_point(aes(x = value, y = bill.length / bill.depth))

summary(lm((bill.length/bill.depth) ~ value, data = df.pheno.adx))

# plot --------------------------------------------------------------------


df.ash.admix.long <- df.ash.admix.long %>% 
  ungroup() %>% 
  #filter(Species == "fortis" | Species == "scandens" | Species == "hybrid") %>% 
  select(sample.ID, Species, nest.Breeding.pop, name, value, First.year.min, Last.year, popGroup = name, prob = value) %>% 
  mutate(nest.Breeding.pop.ext = ifelse(is.na(nest.Breeding.pop), Species, #if breeding pop is NA, assume it is Species ID
                                        ifelse(nest.Breeding.pop == 1, "fortis", #reformat
                                               ifelse(nest.Breeding.pop == 3, "scandens", 
                                                      ifelse(nest.Breeding.pop == 2, "magnirostris", NA)))),
         popGroup.species = ifelse(popGroup == "Cluster2", "magnirostris",
                                   ifelse(popGroup == "Cluster1", "fortis", 
                                          ifelse(popGroup == "Cluster3", "scandens",
                                                 ifelse(popGroup == "Cluster4", "fuliginosa", 
                                                        ifelse(popGroup == "Cluster5", "mag2", NA)))))) %>% 
  select(sample.ID,Species, nest.Breeding.pop.ext, popGroup.species, prob, First.year.min, Last.year)

#####a shortcut to just uses species
#df.fortis_scandens <-  df.fortis_scandens %>% select(-nest.Breeding.pop.ext, nest.Breeding.pop.ext = Species)
#df.fortis_scandens <- df.ash.admix.long %>% select(-Species)
######################

df.ash.admix.wide <- df.ash.admix.long %>% 
  pivot_wider(names_from = popGroup.species, values_from = prob) %>% 
  mutate(ancestry_sca_fortis = scandens - fortis)

df.year.spread <- df.ash.admix.wide %>% 
  filter(nest.Breeding.pop.ext!="hybrid") %>% 
  select(-Species) %>% 
  filter(!is.na(First.year.min)) %>% filter(!is.na(Last.year)) %>% 
  mutate(year = map2(First.year.min, Last.year, `:`)) %>% 
  select(-First.year.min, -Last.year) %>% 
  unnest(cols = c(year))


# classic admixture plot --------------------------------------------------

#build a long form dataframe
cat_ancestry.long <- df.ash.admix.wide %>% 
  select(-First.year.min, -Last.year, -ancestry_sca_fortis) %>% 
  pivot_longer(-c("sample.ID", "nest.Breeding.pop.ext", "Species"), 
               names_to = "popGroup", values_to = "prob") %>% 
  mutate(island = "Daphne", taxa = Species) %>% 
  select(-nest.Breeding.pop.ext, - Species)

#this separates out magnirostris for pretty plottying
#i dont think it was used in the end because of ordering by the dendrogram
cat_ancestry.long <-  cat_ancestry.long%>% 
  mutate(taxa2 = ifelse(taxa!="magnirostris" & taxa!="Unknown", "group1",taxa)) 

cat_ancestry.sum <- cat_ancestry.long %>% 
  group_by(taxa2, popGroup) %>% 
  summarise(mean = mean(prob)) %>% 
  slice_max(mean) %>% 
  select(taxa2, popGroup, max.popGroup = popGroup) %>% 
  left_join(cat_ancestry.long, by = "taxa2")


#various was to rearrange the plot

#colors2 <- c("red","blue","green", "orange", "lightblue")
#
df.ash.admix.wide %>% 
  select(-ancestry_sca_fortis) %>% 
  write_csv("output/admixture/Ash_data/Daphne_admixture_AIMs_Ash_data_K5.csv")
#
#df.4dend <- df.ash.admix.wide %>% select(fortis, magnirostris, scandens, fuliginosa, mag2)
#df.4dend <- sapply(df.4dend, as.numeric )
#
#df.dend <- df.4dend %>% 
#  dist %>% hclust %>% as.dendrogram %>% 
#  set("labels_to_character") %>% color_branches(k=4, col = colors2) %>% set("labels_cex", .05)   
#
#sample.order <- order.dendrogram(df.dend)


#means of ordering 1 where mag is separate and rest sorted by own ancestry
#cat_ancestry.order <- cat_ancestry.sum %>% 
#  filter(popGroup == max.popGroup) %>% 
#  #group_by(taxa2) %>% 
#  arrange(-prob) %>% ungroup() %>% select(sample.ID) %>% 
#  mutate(order = 1:n())

#sorted by dendrogram
#cat_ancestry.order <- cbind(df.ash.admix.wide["sample.ID"], sample.order)
#names(cat_ancestry.order) <- c("sample.ID", "order")

#from relat mat to match figure 1
cat_ancestry.order <- read.csv("output/GEMMA/Daphne_relate_matrix_dendrogram_4k_V3_order.csv")
names(cat_ancestry.order) <- c("sample.ID", "order")


cat_ancestry.sum <- left_join(cat_ancestry.sum, cat_ancestry.order, by = "sample.ID") %>% 
  mutate(popGroup = factor(popGroup, levels = c("mag2","fortis", "magnirostris", "scandens", "fuliginosa", "hybrid")),
         taxa = factor(taxa,  levels = c("mag2","fortis", "magnirostris", "scandens", "fuliginosa", "hybrid")))
palette.colors(palette = "Okabe-Ito")

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'
other.color <-'#D55E00'

cat_ancestry.sum %>% #head(n = 9000) %>% 
  filter(taxa!="Unknown") %>% 
  ggplot(aes(x = fct_reorder(sample.ID, order), y = prob, fill = popGroup)) +
  geom_col(aes(fill = popGroup), alpha = 0.9, size = 0.1) +
  geom_rug(sides = "b", outside = T, aes(color = (taxa)), length = unit(0.060, "npc")) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.0, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_markdown(), #from ggtext for markdown formatting labels
    text = element_text(size = 16)) +
  scale_fill_manual(name = "grp",values = c(other.color, fort.color, magn.color, scan.color,fuli.color), guide = "none") +
  scale_color_manual(name = "Species in the field (bar at bottom)", guide = "none", values = c(fort.color, magn.color, scan.color, fuli.color, hybr.color),
                     labels = c("*G. fortis*", "*G. magnirostris*", "*G. scandens*", "*G. fuliginosa*", "*Geospiza* hybrid")) +
  xlab(NULL)  + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 8)))


ggsave("output/admixture/Ash_data/Daphne_admixture_AIMs_Ash_data_skinny_colors2_k5_pruned.pdf", width = 16, height = 3.5)

# annual plots ------------------------------------------------------------

ann_anc <- function(species, y1, y2, y3, y4){
  
  finch.n <- df.year.spread %>% 
    filter(nest.Breeding.pop.ext == species) %>% 
    group_by(year, nest.Breeding.pop.ext) %>% 
    select(-ancestry_sca_fortis) %>% 
    summarise(n.ancestry = n(), .groups = "drop") %>% 
    select(year, n.ancestry)
  
  finch.geno.freq <- df.year.spread %>% 
    filter(nest.Breeding.pop.ext == species) %>% 
    group_by(year, nest.Breeding.pop.ext) %>% 
    select(-ancestry_sca_fortis) %>% 
    summarise_at(.vars = vars(fortis,magnirostris,scandens,fuliginosa, mag2),
                 .funs = c(mean="mean"))
  
  finch.geno.freq <- left_join(finch.geno.freq, finch.n, by = "year") %>% 
    filter(n.ancestry > 1) 
  
  finch.geno.freq.long <- finch.geno.freq %>% 
    pivot_longer(cols = -c("year", "nest.Breeding.pop.ext", "n.ancestry"), 
                 values_to = "mean.ancestry", names_to = "ancestral.pop")
  
  fort.start <- finch.geno.freq %>% filter(year == 1983)
  
  #title label in italics doesnt work
  #splabel <- paste0("G. ", species)
  #finch.geno.freq.long <- finch.geno.freq.long %>% 
  # mutate(splabel = glue("*{species}*"))
  
  finch.geno.freq.long <- finch.geno.freq.long %>% 
    mutate(ancestral.pop = factor(ancestral.pop, 
                                  levels = c("fortis_mean", "fuliginosa_mean",
                                             "magnirostris_mean", "scandens_mean",
                                             "mag2_mean")))
  
  finch.ancestry <- finch.geno.freq.long %>% 
    filter(year > 1982) %>%
    ggplot(aes(x = year, y = mean.ancestry, color = ancestral.pop)) + 
    geom_hline(yintercept = fort.start$fortis_mean, color = fort.color,  linetype='dotted') +
    geom_hline(yintercept = fort.start$scandens_mean, color = scan.color,  linetype='dotted') +
    geom_hline(yintercept = fort.start$fuliginosa_mean, color = fuli.color,  linetype='dotted') +
    geom_hline(yintercept = fort.start$magnirostris_mean, color = magn.color,  linetype='dotted') +
    geom_hline(yintercept = fort.start$mag2_mean, color = other.color,  linetype='dotted') +
    geom_point(aes(size = n.ancestry)) +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          text = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none") +
    annotate("text",x = 1988, y = y1, label = "G. fortis", fontface = "italic", color = fort.color) +
    annotate("text",x = 1988, y = y2, label = "G. scandens", fontface = "italic", color =  scan.color) +
    annotate("text",x = 1988, y = y3, label = "G. fuliginosa", fontface = "italic", color =  fuli.color) +
    annotate("text",x = 1988, y = y4, label = "G. magnirostris", fontface = "italic", color =  magn.color) +
    scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
    scale_color_manual(values = c(fort.color, fuli.color, magn.color, scan.color, other.color)) +
    labs(y="Ancestry") +
    ylim(0,1)
  
  return(finch.ancestry)

}

#order of labels is fort, scan, ful, mag
p.fort <- ann_anc("fortis", 0.8, 0.125, 0.25, 0.0)
p.fort

p.scan <- ann_anc("scandens", 0.25, 0.85, 0.2, 0.15)
p.scan

p.magn <- ann_anc("magnirostris", -1, -1, -1, 0.85)
p.magn

library(patchwork)
p.fort + p.scan + p.magn
ggsave("output/admixture/Ash_data/annual_change_scandens_and_fortis_all_ancestry_k5.pdf", width = 16, height = 6)


