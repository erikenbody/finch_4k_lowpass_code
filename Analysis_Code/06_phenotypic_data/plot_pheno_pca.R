library(tidyverse)
library(readxl)
library(data.table)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

# load sample list and metadata -------------------------------------------

df.sample.order <- read_table("data/GEMMA/list_of_daphne_vcf_samples.txt", col_names = "meaningful.unique") %>% 
  mutate(order = 1:n())

#load database and 

df.dbm.i <- read_excel("/Users/erikenbody/Dropbox/Shared_Uppsala/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
df.final.proj <- df.dbm.i %>% filter(!is.na(Lowcov_project)) %>% 
  mutate(sex.num = ifelse(is.na(sex.consensus), NA, #convert to numeric
                          ifelse(sex.consensus == "Male", 0, 
                                 ifelse(sex.consensus == "Female", 1, NA)))) 

# merge with VCF sample list ----------------------------------------------
df.final.proj %>% filter(!meaningful.unique %in% df.sample.order$meaningful.unique & Island == "Daphne") %>% nrow()

df.full <- df.sample.order %>% left_join(df.final.proj, by = "meaningful.unique")


df.full2 <- read.csv("output/GEMMA/phenotype_plots/Daphne_phenotype_dataframe_for_plotting.csv")

##drop outliers
df.full <- df.full %>% filter(meaningful.unique!="01Dap21269")
df.full2 <- df.full2 %>% filter(meaningful.unique!="01Dap21269")

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

df.full2 %>% 
  filter(Species == "fortis" | Species == "fuliginosa" | Species == "magnirostris" |
           Species == "scandens" | Species == "hybrid") %>% 
  filter(Island == "Daphne") %>% 
  mutate(Species = factor(Species, levels = c("fortis", "magnirostris", "scandens", "fuliginosa"))) %>% 
  ggplot() + geom_point(aes(x = PC1, y = PC2,
                            color = Species), alpha = 0.8, size = 3) +
  theme_bw() + 
  scale_color_manual(name = "grp",
                     values = c(fort.color, magn.color, scan.color,fuli.color, hybr.color), guide = "none") +
  #scale_color_manual(name = "Species (bar at bottom)", guide = "none", 
  #                   values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')) +
  labs(x = "Beak size (PC1)", y = "Beak shape (PC2)") +
  theme(text = element_text(size = 24))

ggsave("output/GEMMA/phenotype_plots/morphology_PC1_PC2_Daphne_simple.png", width = 8, height = 6)
