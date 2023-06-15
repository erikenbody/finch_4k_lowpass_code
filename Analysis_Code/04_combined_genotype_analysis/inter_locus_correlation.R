library(tidyverse)
library(patchwork)

#setup data. 
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

#drop outliers (two fortis immigrants, one may be maginrostris)
df.all.genos <- df.all.genos %>% 
  filter(meaningful.unique!="01Dap21269" & meaningful.unique!="01Dap21274")

df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")

df.finch.sub <- df.dap.genos %>% filter(Species == "fortis")

tmp.df <- df.finch.sub %>% dplyr::select(all_of(all_loci_linked))
#CramerV(table(tmp.df))
name_fixing <- tibble(locus = names(tmp.df)) # %>% write_csv("output/GEMMA/effect_sizes/tmpnames.csv")

nf2 <- full_join(name_fixing, df_tidy, by = "locus") 
nf2 <- nf2 %>% 
  separate(locus, into = c(NA, NA, "suffix"), remove = F) %>% 
  mutate(suffix = as.numeric(ifelse(grepl("chr", suffix), NA, suffix))) %>% 
  mutate(tidy_locus = ifelse(is.na(tidy_locus) & suffix <10, paste0("G0", suffix), 
                             ifelse(is.na(tidy_locus) & suffix > 10, paste0("G", suffix), tidy_locus))) 

#replace the names of tmp.df with the tidy names 
names(tmp.df) <- nf2$tidy_locus[order(match(nf2$locus, names(tmp.df)))]

df_numeric <- tmp.df

df_numeric[tmp.df == "AA" ] <- 0
df_numeric[tmp.df == "AB" ] <- 1
df_numeric[tmp.df == "BB" ] <- 2
df_numeric[] <- lapply(df_numeric, as.numeric)

calculate_r2 <- function(x, y) {
  correlation <- cor(x, y, use = "pairwise.complete.obs")  
  return(correlation^2)
}


ld_matrix <- matrix(0, ncol = ncol(df_numeric), nrow = ncol(df_numeric))
for (i in 1:(ncol(df_numeric)-1)) {
  for (j in (i+1):ncol(df_numeric)) {
    ld_matrix[i, j] <- calculate_r2(df_numeric[,i], df_numeric[,j])
    ld_matrix[j, i] <- ld_matrix[i, j]  # LD is symmetric
  }
}


library(pheatmap)
# Get column names from original dataframe
colnames(ld_matrix) <- colnames(df_numeric)
rownames(ld_matrix) <- colnames(df_numeric)

all_loci_linked.order <- c("G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09","G29", "G10", "G11", "G12", "G13",
                           "G14", "G15", "G16", "G17", "G18", "G19","G30", "G20", "G21", "G22", "G23", "G24", "G25", "G26",
                           "G27", "G28")

ld_matrix_ordered <- ld_matrix[all_loci_linked.order, all_loci_linked.order]
diag(ld_matrix_ordered) <- 1
ld_matrix_ordered_rev <- ld_matrix_ordered[rev(rownames(ld_matrix_ordered)), ]

# Create heatmap
pheatmap(ld_matrix_ordered_rev, 
         main = "Linkage Disequilibrium (R2) Heatmap", 
         color = colorRampPalette(c("white", "red"))(25),  # gradient of colors from white to red
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_names_col = TRUE)

library(reshape2)


ld_matrix_df <- as.data.frame(ld_matrix_ordered_rev)
ld_matrix_df$locus <- rownames(ld_matrix_ordered_rev)

# Melt the data frame to a long format
ld_matrix_melt <- melt(ld_matrix_df, id.vars="locus")
colnames(ld_matrix_melt) <- c("locus1", "locus2", "R2")

ld_matrix_melt$locus1 <- factor(ld_matrix_melt$locus1, levels = all_loci_linked.order)
ld_matrix_melt$locus2 <- factor(ld_matrix_melt$locus2, levels = all_loci_linked.order)  # Reverse the order for the y-axis


# Create the heatmap
ggplot(ld_matrix_melt, aes(x=locus1, y=locus2, fill=R2)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue",
                      guide = guide_colourbar(barwidth = 5, barheight = 30)) +  
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 22),
        axis.text.y=element_text(size = 22),
        legend.text=element_text(size=22),   # Adjust the size of the legend text
        legend.title=element_text(size=22)) +
  labs(x="", y="", fill="r^2") +
  geom_text(aes(label=round(R2, 2)), size=6) 
ggsave("output/GEMMA/genotype_plots/r2_genotype_all_loci.png", width = 20, height = 20)

