library(tidyverse)
library(RColorBrewer)
library(heatmap3)
library(dendextend)

df.relate.matrix <- read.table("data/GEMMA/Daphne_3957/autosomes_noNA.relate.cXX.txt")
mat.relate.matrix <- as.matrix(df.relate.matrix)

sample.names.in <- read_table("data/GEMMA/list_of_daphne_vcf_samples.txt", col_names = "meaningful.unique") 
colnames(df.relate.matrix) <- sample.names.in$meaningful.unique
rownames(df.relate.matrix) <- sample.names.in$meaningful.unique

df.dbm.i <- read_excel("/Users/erikenbody/Dropbox/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
df.final.proj <- df.dbm.i %>% filter(!is.na(Lowcov_project)) 

sample.names <- left_join(sample.names.in, df.final.proj, by = "meaningful.unique")

sample.names$colours <- colorRampPalette(brewer.pal(7, "Pastel1"))(7)[as.factor(sample.names$Species)]
sample.names$colours.species <- colorRampPalette(brewer.pal(7, "Pastel1"))(7)[as.factor(sample.names$Species)]

# dendrogram --------------------------------------------------------------

dend <- df.relate.matrix %>% dist %>% hclust %>% as.dendrogram %>%  set("labels_to_character") %>% color_branches(k=3)

df.dend.labs <- data.frame(meaningful.unique = labels(dend))
df.dend.labs <- left_join(df.dend.labs, sample.names, by = c("meaningful.unique" = "meaningful.unique"))
df.dend.labs %>% 
  mutate(order = 1:n()) %>% 
  select(meaningful.unique, order) %>% 
  write_csv("output/GEMMA/Daphne_relate_matrix_dendrogram_4k_V3_order.csv")

labels_colors(dend) <- df.dend.labs$colours
#new colors
dend <- dend %>% set("branches_k_color", value = c(scan.color,magn.color, fort.color), k = 3)

dend <- dend %>%  set("labels_cex", .05)   

pdf("output/GEMMA/Daphne_relate_matrix_dendrogram_4k_V2.pdf", width = 18, height = 10)
plot(dend)
dev.off()

pdf("output/GEMMA/Daphne_relate_matrix_dendrogram_4k_V3.pdf", width = 16, height = 4)
labels_colors(dend) <- "white"
plot(dend)
dev.off()


clusters <- cutree(dend, k=3)

cluster1 <- names(clusters[clusters == "1"]) #fortis 01
cluster2 <- names(clusters[clusters == "2"]) #scandens 03
cluster3 <- names(clusters[clusters == "3"]) #mag 02

grep("01Dap", cluster1, value = T, invert = T)
grep("03Dap", cluster2, value = T, invert = T)
grep("02Dap", cluster3, value = T, invert = T)

#cluster 1 is mostly fortis, cluster 2 is mostly scandens, cluster 3 is mag
write_csv(as.data.frame(cluster1), "output/GEMMA/Daphne_relatedness_cluster_1.csv", col_names = F)
write_csv(as.data.frame(cluster2), "output/GEMMA/Daphne_relatedness_cluster_2.csv", col_names = F)
write_csv(as.data.frame(cluster3), "output/GEMMA/Daphne_relatedness_cluster_3.csv", col_names = F)

