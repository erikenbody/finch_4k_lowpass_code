library(data.table)
library("tidyverse")
library(ggrepel)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(cluster)

# using all samples -------------------------------------------------------


df.dbm.i2 <- read_excel("/Users/erikenbody/Dropbox/Shared_Uppsala/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
df.dbm.2 <- df.dbm.i2 %>% filter(!is.na(Lowcov_project) | !is.na(High_coverage_project))
#these have wrong names because they were generated accidentally originally.
#whole genome all sample
full.pca <- fread("data/plink_pca/all_highcov_lowcov/chr5_region_Daphne_3species.eigenvec", header = TRUE)
full.val <- fread("data/plink_pca/all_highcov_lowcov/chr5_region_Daphne_3species.eigenval", header = FALSE)
full.val$prop <- (full.val$V1 / (sum(full.val$V1))) * 100

i.full.pca <- full.pca %>% left_join(df.dbm.2, by = c("#IID" = "meaningful.unique")) %>% 
  dplyr::select(`#IID`,TubeNum, High_coverage_project,Island,Genus, Page.01, Slot.02, NUM.07, AltNum.08, mother, social.father, EP.father, Symbol, UU_notes, `DNA.dilution.plate…`, `Pool (2020)`, First.year.min, Last.year, genetic.sex, Genotype, Category, Species, starts_with("msatt"), starts_with("PC"), starts_with("bill"), meaningful.unique = `#IID`) %>% 
  mutate(species.ext = ifelse(Species == "hybrid" & !grepl("f", Genotype), "fortisxscandens",
                              ifelse(Species == "hybrid" & !grepl("S", Genotype), "fortisxfuliginosa", 
                                     ifelse(Species == "hybrid", "trihybrid", Species))),
         geospize.sp = ifelse(Genus =="Geospiza", Species, NA))

i.full.pca %>% filter(PC1 <0 & PC2 < -0.02 & Species!="magnirostris")

p1 <- ggplot(i.full.pca) + geom_point(aes(x = PC1, y = PC2, color = geospize.sp), alpha = 0.8) + theme_bw() +
  theme(text = element_text(size = 14)) +
  #geom_text_repel(segment.size = 0.05, size = 9, aes(label = label2, color = species.grouping)) +
  xlab(paste0("PC1", ": ", round(full.val[1,2],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(full.val[2,2],1),"% variance"))

p1

i.full.pca %>% filter(Species == "magnirostris") %>% 
  ggplot() + geom_point(aes(x = PC1, y = PC2, color = Island), alpha = 0.8, size = 4) + theme_bw() +
  theme(text = element_text(size = 14)) +
  #geom_text_repel(segment.size = 0.05, size = 9, aes(label = label2, color = species.grouping)) +
  xlab(paste0("PC1", ": ", round(full.val[1,2],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(full.val[2,2],1),"% variance")) +
  scale_colour_brewer(palette = "Set2")
ggsave("output/plink_pca/magnirostris_analysis/magnirostris_from_whole_sample_PCA.png", width = 12, height = 10)

i.full.pca %>% filter(Species == "magnirostris" | Species == "fortis") %>% 
  ggplot() + geom_point(aes(x = PC1, y = PC2, color = Island, shape = Species), alpha = 0.8, size = 4) + theme_bw() +
  theme(text = element_text(size = 14)) +
  #geom_text_repel(segment.size = 0.05, size = 9, aes(label = label2, color = species.grouping)) +
  xlab(paste0("PC1", ": ", round(full.val[1,2],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(full.val[2,2],1),"% variance"))

i.full.pca %>% mutate(num.num = as.numeric(TubeNum)) %>% 
  filter(num.num %in% c(6682,6684,9283,9389))

i.full.pca %>% filter(Species == "magnirostris") %>% 
  mutate(num.num = as.numeric(TubeNum)) %>% 
  mutate(check.these = ifelse(num.num %in% c(6682,6684,9283,9389), "yes","no")) %>% 
  ggplot() + geom_point(aes(x = PC1, y = PC2, color = check.these))

# clustering --------------------------------------------------------------

mag.pca <- i.full.pca %>% filter(Species == "magnirostris") 
set.seed(42)
k <- kmeans(mag.pca[,c("PC1","PC2")], 2, nstart = 25, iter.max = 100)
sil <- silhouette(k$cluster, dist(mag.pca[,c("PC1","PC2")]))

##Assign the cluster value and the confidence value to the main dataframe
mag.pca$cluster <- paste("cluster",k$cluster, sep = "_")
mag.pca$sil <- sil[,3]

##assign a threshold. In this case, values < 0.5 clearly cluster intermediately between homozygous and heterozygous individuals
mag.pca <- mag.pca %>% mutate(cluster.pca = ifelse(sil < .1, NA, cluster))


mag.pca$cluster <- paste("cluster",k$cluster, sep = "_")

mag.pca %>% 
  ggplot() + geom_point(aes(x = PC1, y = PC2, color = cluster.pca), alpha = 0.8, size = 4) + theme_bw() +
  theme(text = element_text(size = 14))

mag.pca <- mag.pca %>% mutate(grouping = ifelse(PC2 < -0.03 & PC1 < -0.007, "A","B")) 
ggplot(mag.pca) + geom_point(aes(x = PC1, y = PC2, color = grouping), alpha = 0.8, size = 4) + theme_bw() +
  theme(text = element_text(size = 14))

mag.pca <- mag.pca %>% mutate(Daphne.grouping = ifelse(grouping == "A", "Daphne_A",
                                            ifelse(grouping == "B" & Island == "Daphne", "Daphne_B", Island)))

ggplot(mag.pca) + 
  geom_point(aes(x = PC1, y = PC2, color = Daphne.grouping), alpha = 0.8, size = 4) + theme_bw() +
  theme(text = element_text(size = 14))
ggsave("output/plink_pca/magnirostris_analysis/Daphne_magnirostris_groupings_PCA.png", width = 12, height = 10)

mag.pca.out <- mag.pca %>% 
  dplyr::select(Daphne.grouping, meaningful.unique, TubeNum, Page.01, Slot.02, NUM.07, AltNum.08, 
         Species, Island, bill.depth, bill.width, bill.length, 
         First.year.min, Last.year, genetic.sex, mother, social.father, EP.father, PC1, PC2)
#dir.create("output/plink_pca/magnirostris_analysis")
write.csv(mag.pca.out, "output/plink_pca/magnirostris_analysis/magnirostris_Daphne_groupings.csv", na = "")
mag.pca.out %>% filter(Island == "Daphne") %>% 
  dplyr::select(meaningful.unique, Daphne.grouping) %>% 
  write_tsv("output/plink_pca/magnirostris_analysis/Daphne_magnirostris_groupings_for_pixy.txt", col_names = F)

mag.pca.out %>% filter(PC1 > 0) 
mag.pca.out %>% filter(PC1 < -0.005 & PC2 > -0.02)


mag.pca.out %>% filter(Island == "Daphne") %>% 
  ggplot() +
  geom_jitter(aes(x = Daphne.grouping, y = bill.depth + bill.width + bill.length), width = 0.1, alpha = 0.5) +
  geom_boxplot(aes(x = Daphne.grouping, y = bill.depth + bill.width + bill.length), fill = NA, outlier.shape = NA) +
    theme_bw() + xlab("magnirostris groupings on Daphne")
ggsave("output/plink_pca/magnirostris_analysis/Daphne_magnirostris_groupings_phenotype.png", width = 8, height = 6)

mag.pca.dap <- mag.pca.out %>% filter(Island == "Daphne")
t.test(mag.pca.dap$bill.depth ~ mag.pca.dap$Daphne.grouping)
mag.pca.dap %>% 
  group_by(Daphne.grouping) %>% 
  summarise(meansize = mean((bill.depth + bill.width + bill.length), na.rm = T))

# -------------------------------------------------------------------------

# magnirostris ------------------------------------------------------------

#only use high coverage ref panel
mag_samps <- mag.pca %>% 
  filter(Island == "Daphne" &Species == "magnirostris" | 
           Island!="Daphne" & Species == "magnirostris" & !is.na(High_coverage_project)) %>% 
  filter(meaningful.unique!="02Dap9389") #remove a sample that is likely some kind of immigrant hybrid

mag_samps %>% 
  filter(Species == "magnirostris") %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = Daphne.grouping), alpha = 0.8, size = 4) + 
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.title = element_blank()) +
  #geom_text_repel(segment.size = 0.05, size = 9, aes(label = label2, color = species.grouping)) +
  xlab(paste0("PC1", ": ", round(full.val[1,2],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(full.val[2,2],1),"% variance")) +
  scale_colour_brewer(palette = "Set2") 
ggsave("output/plink_pca/magnirostris_from_whole_sample_PCA.png", width = 12, height = 10)

