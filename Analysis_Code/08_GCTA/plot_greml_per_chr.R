library(tidyverse)

chr_order <- c("chr1", "chr1A" , "chr2", "chr3", "chr4","chr4A", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11","chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
               "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chrZ",  "chrunknown")

df.pchr <- read.table("data/gcta/compare_all_chr_greml_VGVP.txt")
df.pchr2 <- read.table("data/gcta/compare_all_chr_greml_VGVP_incluZ.txt")
plot(df.pchr$V2 ~ df.pchr2$V2[1:30])
df.chridx <- read.table("data/gcta/full_ordered_chr.list")
df.chridx<-df.chridx %>% mutate(V1 = ifelse(V1 == "chrLGE22", "chr29", V1))
df.chridx_noZ <- df.chridx %>% filter(V1!="chrZ")

df.pchr2$chr <- factor(df.chridx$V1, levels = chr_order)

ggplot(data = df.pchr2, aes(x = chr, y = V2)) +
  geom_point() +
  geom_line(group = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "Vg/Vp")
ggsave("output/gcta/greml_per_chr.png", width = 10, height = 3)
