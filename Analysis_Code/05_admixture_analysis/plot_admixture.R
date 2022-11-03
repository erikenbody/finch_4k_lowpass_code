library(tidyverse)
library(patchwork)
library(readxl)
library(dendextend)
library(ggtext)

#df.ash.admix <- read.table("data/admixture/Daphne_FortisFuliginosaScandens_K3_AdmixProportions.txt", header = T)
df.ash.admix_og <- read.table("data/admixture/Ash_data/with_hybrids/projAdmix.4.Q", header = F)
names(df.ash.admix_og) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4")
sample.names <- read.table("data/admixture/Ash_data/with_hybrids/sample_order.txt", header = F)
df.ash.admix_og$sample.ID <- sample.names$V1

df.ash.admix <- read.table("data/admixture/Ash_data/admix_proportions_PRUNED/K4_proportions.txt", header = T)
df.ash.admix$sample.ID <- df.ash.admix$meaningful.unique ; df.ash.admix$meaningful.unique <- NULL


V2.mdb <- read_excel("/Users/erikenbody/Dropbox/Shared_Uppsala/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)

df.final.proj <- V2.mdb %>% filter(!is.na(Lowcov_project) | !is.na(High_coverage_project))
df.final.proj %>% filter(Genus == "Geospiza" & meaningful.unique!="02Gen1389") %>%
  select(meaningful.unique, Species,Island, Genotype, nest.Breeding.pop, contains("bill")) %>%
  distinct(meaningful.unique, .keep_all = T) %>%
  write_csv("data/snmf/geospiza_population_info.csv")

pop_info <- df.final.proj %>% filter(Genus == "Geospiza" & meaningful.unique!="02Gen1389") %>%
  #select(meaningful.unique) %>%
  select(meaningful.unique, Species,Island, Genotype, nest.Breeding.pop, contains("bill"), First.year.min, Last.year) %>%
  distinct(meaningful.unique, .keep_all = T)

df.ash.admix.long <- df.ash.admix %>% 
  pivot_longer(cols = -sample.ID)

df.ash.admix.long <- left_join(df.ash.admix.long, pop_info, by = c("sample.ID" = "meaningful.unique"))

df.ash.admix.long %>% 
  group_by(Species, name) %>% 
  summarise(mean.val = mean(value)) %>% 
  print(n=100)


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
                                                 ifelse(popGroup == "Cluster4", "fuliginosa", NA))))) %>% 
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

cat_ancestry.long <- df.ash.admix.wide %>% 
  select(-First.year.min, -Last.year, -ancestry_sca_fortis) %>% 
  pivot_longer(-c("sample.ID", "nest.Breeding.pop.ext", "Species"), 
               names_to = "popGroup", values_to = "prob") %>% 
  mutate(island = "Daphne", taxa = Species) %>% 
  select(-nest.Breeding.pop.ext, - Species)
  
#this separates out magnirostris for pretty plottying
cat_ancestry.long <-  cat_ancestry.long%>% 
  mutate(taxa2 = ifelse(taxa!="magnirostris" & taxa!="Unknown", "group1",taxa)) 

cat_ancestry.sum <- cat_ancestry.long %>% 
  group_by(taxa2, popGroup) %>% 
  summarise(mean = mean(prob)) %>% 
  slice_max(mean) %>% 
  select(taxa2, popGroup, max.popGroup = popGroup) %>% 
  left_join(cat_ancestry.long, by = "taxa2")

colors2 <- c("red","blue","green", "orange")


df.4dend <- df.ash.admix.wide %>% select(fortis, magnirostris, scandens, fuliginosa)
df.4dend <- sapply(df.4dend, as.numeric )

df.dend <- df.4dend %>% 
  dist %>% hclust %>% as.dendrogram %>% 
  set("labels_to_character") %>% color_branches(k=4, col = colors2) %>% set("labels_cex", .05)   

sample.order <- order.dendrogram(df.dend)


#means of ordering 1 where mag is separate and rest sorted by own ancestry
cat_ancestry.order <- cat_ancestry.sum %>% 
  filter(popGroup == max.popGroup) %>% 
  #group_by(taxa2) %>% 
  arrange(-prob) %>% ungroup() %>% select(sample.ID) %>% 
  mutate(order = 1:n())

#sorted by dendrogram
#cat_ancestry.order <- cbind(df.ash.admix.wide["sample.ID"], sample.order)
#names(cat_ancestry.order) <- c("sample.ID", "order")

#from relat mat 
cat_ancestry.order <- read.csv("output/GEMMA/Daphne_relate_matrix_dendrogram_4k_V3_order.csv")
names(cat_ancestry.order) <- c("sample.ID", "order")


cat_ancestry.sum <- left_join(cat_ancestry.sum, cat_ancestry.order, by = "sample.ID") %>% 
  mutate(popGroup = factor(popGroup, levels = c("fortis", "magnirostris", "scandens", "fuliginosa", "hybrid")),
         taxa = factor(taxa,  levels = c("fortis", "magnirostris", "scandens", "fuliginosa", "hybrid")))
palette.colors(palette = "Okabe-Ito")

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

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
  scale_fill_manual(name = "grp",values = c(fort.color, magn.color, scan.color,fuli.color), guide = "none") +
  scale_color_manual(name = "Species (bar at bottom)", guide = "none", values = c(fort.color, magn.color, scan.color, fuli.color, hybr.color),
                     labels = c("*G. fortis*", "*G. magnirostris*", "*G. scandens*", "*G. fuliginosa*", "*Geospiza* hybrid")) +
  xlab(NULL)  + 
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 8)))

ggsave("output/admixture/Ash_data/Daphne_admixture_AIMs_Ash_data_skinny_colors2_v2.pdf", width = 16, height = 3.5)


# plot just scandens and fortis -------------------------------------------

###explor#
df.year.spread %>% filter(nest.Breeding.pop.ext == "fortis") %>%
  group_by(year) %>% 
  summarise(mean_anc = mean(scandens, na.rm = T), sd_anc = sd(scandens), n = n()) %>% print(n = 50)
#ggplot() + geom_point(aes(x = year, y = mean_anc))

df.year.spread %>% filter(nest.Breeding.pop.ext == "scandens") %>%
  group_by(year) %>% 
  summarise(mean_anc = mean(fortis, na.rm = T), sd_anc = sd(fortis), n = n()) %>% print(n = 50)
###
#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

adult.geno.freq <- df.year.spread %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  # mutate(ancestry.01 = range01(ancestry_sca_fortis)) %>% 
  summarise(mean.ancestry = mean(ancestry_sca_fortis, na.rm = TRUE),
            sd.ancestry_sca_fortis = sd(ancestry_sca_fortis),
            n.ancestry = n()) %>%
  filter(n.ancestry > 1) %>% 
  mutate(se.ancestry = sd.ancestry_sca_fortis / sqrt(n.ancestry),
         lower.ci.ancestry = mean.ancestry - qt(1 - (0.05 / 2), n.ancestry - 1) * se.ancestry,
         upper.ci.ancestry = mean.ancestry + qt(1 - (0.05 / 2), n.ancestry - 1) * se.ancestry)

df.scan <- adult.geno.freq %>% filter(nest.Breeding.pop.ext == "scandens" & year == 1983)
df.fort <- adult.geno.freq %>% filter(nest.Breeding.pop.ext == "fortis" & year == 1983)

adult.geno.freq %>% 
  filter(nest.Breeding.pop.ext == "scandens" | nest.Breeding.pop.ext == "fortis") %>% 
  filter(year > 1982) %>%
  ggplot(aes(x = year, y = mean.ancestry, color = nest.Breeding.pop.ext)) + 
  geom_hline(yintercept = df.scan$mean.ancestry, color = "#67a9cf",  linetype='dotted') +
  geom_hline(yintercept = df.fort$mean.ancestry, color = "#ef8a62",  linetype='dotted') +
  geom_hline(yintercept = 0.5, color = "black",  linetype='dotted') +
  geom_hline(yintercept = -0.5, color = "black",  linetype='dotted') +
  geom_point(aes(size = n.ancestry)) +
  geom_line() +
  #geom_linerange(aes(ymin = lower.ci.ancestry, ymax = upper.ci.ancestry), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text",x = 1988, y = -0.7, label = "G. fortis", fontface = "italic", color = "#ef8a62") +
  annotate("text",x = 1988, y = 0.7, label = "G. scandens", fontface = "italic", color =  "#67a9cf") +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1983, 2012, by = 1)) +
  scale_color_manual(values = c("#ef8a62", "#67a9cf")) +
  labs(y=expression(paste(italic("G. scandens")," ancestry - ", italic("G.fortis"), " ancestry"))) +
  ylim(-1,1)

#ggsave("output/snmf/annual_change_scandens_fortis_ancestry_V3.pdf", width = 8, height = 6)
ggsave("output/admixture/Ash_data/annual_change_scandens_fortis_ancestry_V4.pdf", width = 8, height = 6)


# fortis all ancestry ---------------------------------------------------

fortis.n <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "fortis") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise(n.ancestry = n(), .groups = "drop") %>% 
  select(year, n.ancestry)

fortis.geno.freq <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "fortis") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise_at(.vars = vars(fortis,magnirostris,scandens,fuliginosa),
               .funs = c(mean="mean"))
fortis.geno.freq <- left_join(fortis.geno.freq, fortis.n, by = "year") %>% 
  filter(n.ancestry > 1) 

fortis.geno.freq.long <- fortis.geno.freq %>% 
  pivot_longer(cols = -c("year", "nest.Breeding.pop.ext", "n.ancestry"), 
               values_to = "mean.ancestry", names_to = "ancestral.pop")

fort.start <- fortis.geno.freq %>% filter(year == 1983)

fortis.ancestry <- fortis.geno.freq.long %>% 
  filter(year > 1982) %>%
  ggplot(aes(x = year, y = mean.ancestry, color = ancestral.pop)) + 
  geom_hline(yintercept = fort.start$fortis_mean, color = fort.color,  linetype='dotted') +
  geom_hline(yintercept = fort.start$scandens_mean, color = scan.color,  linetype='dotted') +
  geom_hline(yintercept = fort.start$fuliginosa_mean, color = fuli.color,  linetype='dotted') +
  geom_hline(yintercept = fort.start$magnirostris_mean, color = hybr.color,  linetype='dotted') +
  #geom_hline(yintercept = 0.5, color = "black",  linetype='dotted') +
  #geom_hline(yintercept = -0.5, color = "black",  linetype='dotted') +
  geom_point(aes(size = n.ancestry)) +
  geom_line() +
  #geom_linerange(aes(ymin = lower.ci.ancestry, ymax = upper.ci.ancestry), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text",x = 1988, y = 0.8, label = "G. fortis", fontface = "italic", color = fort.color) +
  annotate("text",x = 1988, y = 0.0, label = "G. scandens", fontface = "italic", color =  scan.color) +
  annotate("text",x = 1988, y = 0.2, label = "G. fuliginosa", fontface = "italic", color =  fuli.color) +
  annotate("text",x = 1988, y = 0.1, label = "G. magnirostris", fontface = "italic", color =  magn.color) +
  #scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1983, 2012, by = 1)) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  scale_color_manual(values = c(fort.color,fuli.color,magn.color,scan.color)) +
  #labs(y=expression(paste("Ancestry in ", italic("G. fortis")))) +
  labs(y="Ancestry") +
  ylim(0,1) +
  ggtitle(expression(paste(italic("G. fortis"))))


#fortis.ancestry
write_csv(fortis.geno.freq.long, "output/admixture/fortis_admixture_longform_annual.csv")

# scandens all ancestry ---------------------------------------------------

scandens.n <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "scandens") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise(n.ancestry = n(), .groups = "drop") %>% 
  select(year, n.ancestry)

scandens.geno.freq <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "scandens") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise_at(.vars = vars(fortis,magnirostris,scandens,fuliginosa),
               .funs = c(mean="mean"))
scandens.geno.freq <- left_join(scandens.geno.freq, scandens.n, by = "year") %>% 
  filter(n.ancestry > 1) 

scandens.geno.freq.long <- scandens.geno.freq %>% 
  pivot_longer(cols = -c("year", "nest.Breeding.pop.ext", "n.ancestry"), 
               values_to = "mean.ancestry", names_to = "ancestral.pop")

scan.start <- scandens.geno.freq %>% filter(year == 1983)

scandens.ancestry <- scandens.geno.freq.long %>% 
  filter(year > 1982) %>%
  ggplot(aes(x = year, y = mean.ancestry, color = ancestral.pop)) + 
  geom_hline(yintercept = scan.start$fortis_mean, color = fort.color,  linetype='dotted') +
  geom_hline(yintercept = scan.start$scandens_mean, color = scan.color,  linetype='dotted') +
  geom_hline(yintercept = scan.start$fuliginosa_mean, color = fuli.color,  linetype='dotted') +
  geom_hline(yintercept = scan.start$magnirostris_mean, color = hybr.color,  linetype='dotted') +
  #geom_hline(yintercept = 0.5, color = "black",  linetype='dotted') +
  #geom_hline(yintercept = -0.5, color = "black",  linetype='dotted') +
  geom_point(aes(size = n.ancestry)) +
  geom_line() +
  #geom_linerange(aes(ymin = lower.ci.ancestry, ymax = upper.ci.ancestry), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text",x = 1985, y = 0.19, label = "G. fortis", fontface = "italic", color = fort.color) +
  annotate("text",x = 1985, y = 0.87, label = "G. scandens", fontface = "italic", color =  scan.color) +
  annotate("text",x = 1985, y = 0.16, label = "G. fuliginosa", fontface = "italic", color =  fuli.color) +
  annotate("text",x = 1985, y = 0.13, label = "G. magnirostris", fontface = "italic", color =  magn.color) +
  #scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1983, 2012, by = 1)) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  scale_color_manual(values = c(fort.color,fuli.color,magn.color,scan.color)) +
  labs(y=expression(paste("Ancestry in ", italic("G. scandens")))) +
  ylim(0,1) +
  ggtitle(expression(paste(italic("G. scandens"))))+ ylab(NULL)
#scandens.ancestry

write_csv(scandens.geno.freq.long, "output/admixture/scandens_admixture_longform_annual.csv")


# magnirostris ancestry ---------------------------------------------------

magnirostris.n <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "magnirostris") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise(n.ancestry = n(), .groups = "drop") %>% 
  select(year, n.ancestry)

magnirostris.geno.freq <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "magnirostris") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise_at(.vars = vars(fortis,scandens,magnirostris,fuliginosa),
               .funs = c(mean="mean"))
magnirostris.geno.freq <- left_join(magnirostris.geno.freq, magnirostris.n, by = "year") %>% 
  filter(n.ancestry > 1) 

magnirostris.geno.freq.long <- magnirostris.geno.freq %>% 
  pivot_longer(cols = -c("year", "nest.Breeding.pop.ext", "n.ancestry"), 
               values_to = "mean.ancestry", names_to = "ancestral.pop")

magn.start <- magnirostris.geno.freq %>% filter(year == 1984)

magnirostris.ancestry <- magnirostris.geno.freq.long %>% 
  filter(year > 1982) %>%
  ggplot(aes(x = year, y = mean.ancestry, color = ancestral.pop)) + 
  geom_hline(yintercept = magn.start$fortis_mean, color = fort.color,  linetype='dotted') +
  geom_hline(yintercept = magn.start$scandens_mean, color = scan.color,  linetype='dotted') +
  geom_hline(yintercept = magn.start$fuliginosa_mean, color = fuli.color,  linetype='dotted') +
  geom_hline(yintercept = magn.start$magnirostris_mean, color = hybr.color,  linetype='dotted') +
  #geom_hline(yintercept = 0.5, color = "black",  linetype='dotted') +
  #geom_hline(yintercept = -0.5, color = "black",  linetype='dotted') +
  geom_point(aes(size = n.ancestry)) +
  geom_line() +
  #geom_linerange(aes(ymin = lower.ci.ancestry, ymax = upper.ci.ancestry), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text",x = 1984.5, y = 0.15, label = "G. fortis", fontface = "italic", color = fort.color) +
  annotate("text",x = 1984.5, y = 0.03, label = "G. scandens", fontface = "italic", color =  scan.color) +
  annotate("text",x = 1984.5, y = 0.05, label = "G. fuliginosa", fontface = "italic", color =  fuli.color) +
  annotate("text",x = 1985, y = 0.95, label = "G. magnirostris", fontface = "italic", color =  magn.color) +
  #scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1983, 2012, by = 1)) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  scale_color_manual(values = c(fort.color, fuli.color, magn.color, scan.color)) +
  labs(y=expression(paste("Ancestry in ", italic("G. magnirostris")))) +
  ylim(0,1) +
  ggtitle(expression(paste(italic("G. magnirostris"))))+ ylab(NULL)

magnirostris.ancestry

# fuliginosa ancestry ---------------------------------------------------

fuliginosa.n <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "fuliginosa") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise(n.ancestry = n(), .groups = "drop") %>% 
  select(year, n.ancestry)

fuliginosa.geno.freq <- df.year.spread %>% 
  filter(nest.Breeding.pop.ext == "fuliginosa") %>% 
  group_by(year, nest.Breeding.pop.ext) %>% 
  select(-ancestry_sca_fortis) %>% 
  summarise_at(.vars = vars(fortis,scandens,magnirostris,fuliginosa),
               .funs = c(mean="mean"))
fuliginosa.geno.freq <- left_join(fuliginosa.geno.freq, fuliginosa.n, by = "year") %>% 
  filter(n.ancestry > 1) 

fuliginosa.geno.freq.long <- fuliginosa.geno.freq %>% 
  pivot_longer(cols = -c("year", "nest.Breeding.pop.ext", "n.ancestry"), 
               values_to = "mean.ancestry", names_to = "ancestral.pop")

fuli.start <- fuliginosa.geno.freq %>% filter(year == 1983)

fuliginosa.ancestry <- fuliginosa.geno.freq.long %>% 
  filter(year > 1982) %>%
  ggplot(aes(x = year, y = mean.ancestry, color = ancestral.pop)) + 
  geom_hline(yintercept = fuli.start$fortis_mean, color = fort.color,  linetype='dotted') +
  geom_hline(yintercept = fuli.start$scandens_mean, color = scan.color,  linetype='dotted') +
  geom_hline(yintercept = fuli.start$fuliginosa_mean, color = fuli.color,  linetype='dotted') +
  geom_hline(yintercept = fuli.start$magnirostris_mean, color = hybr.color,  linetype='dotted') +
  #geom_hline(yintercept = 0.5, color = "black",  linetype='dotted') +
  #geom_hline(yintercept = -0.5, color = "black",  linetype='dotted') +
  geom_point(aes(size = n.ancestry)) +
  geom_line() +
  #geom_linerange(aes(ymin = lower.ci.ancestry, ymax = upper.ci.ancestry), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text",x = 1984, y = 0.15, label = "G. fortis", fontface = "italic", color = fort.color) +
  annotate("text",x = 1984, y = 0.03, label = "G. scandens", fontface = "italic", color =  scan.color) +
  annotate("text",x = 1984, y = 0.05, label = "G. fuliginosa", fontface = "italic", color =  fuli.color) +
  annotate("text",x = 1984, y = 0.95, label = "G. magnirostris", fontface = "italic", color =  magn.color) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1983, 2012, by = 1)) +
  scale_color_manual(values = c(fort.color, fuli.color, magn.color, scan.color)) +
  labs(y=expression(paste("Ancestry in ", italic("G. fuliginosa")))) +
  ylim(0,1) +
  ggtitle(expression(paste(italic("G. fuliginosa")))) + ylab(NULL)

fuliginosa.ancestry


# output plots ------------------------------------------------------------

fortis.ancestry + scandens.ancestry
ggsave("output/admixture/Ash_data/annual_change_scandens_and_fortis_all_ancestry_V2.pdf", width = 14, height = 8.5)

fortis.ancestry + scandens.ancestry + magnirostris.ancestry
ggsave("output/admixture/Ash_data/annual_change_scandens_and_fortis_all_ancestry_V2.pdf", width = 16, height = 6)


# ancestry summaries ------------------------------------------------------

fortis.geno.freq.long %>% filter(year == 1983 | year == 2012)

scandens.geno.freq.long %>% filter(year == 1983 | year == 2012)

magnirostris.geno.freq.long %>% filter(year == 1984 | year == 2012)
