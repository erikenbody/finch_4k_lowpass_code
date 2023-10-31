library(tidyverse)

df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

df.mag <- df.all.genos %>% 
  filter(Species == "magnirostris" & Island == "Daphne")

df.mag <- df.mag %>% 
  select(meaningful.unique, bill.PC1.sp, bill.depth, bill.width, bill.length, anc_mag, anc_mag2, magnirostris_grouping) %>% 
  mutate(bill.size = bill.depth + bill.width + bill.length)

#eye balled grouping from PCA analysis
#df.mag %>% 
#  ggplot() + 
#  geom_point(aes(x = anc_mag, y = anc_mag2, color = magnirostris_grouping))

#set a threshold: if > 50% ancestry from mag2, consider it mag2

df.mag <- df.mag %>% 
  mutate(mag_admx_group = ifelse(anc_mag2 > 0.5, "Group_B", "Group_A")) %>% 
  filter(!is.na(mag_admx_group))

df.mag %>% 
  ggplot() + 
  geom_point(aes(x = anc_mag, y = anc_mag2, color = mag_admx_group))

ggplot(df.mag, aes(x = mag_admx_group, y = bill.depth + bill.width + bill.length)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  xlab("Daphne magnirostris grouping")

ggsave("output/admixture/Ash_data/mag_two_clusters_bill_size.png", width = 6, height = 5)

test_result <- t.test(df.mag$bill.size ~ df.mag$mag_admx_group)
47.65682-45.23842 

2.4184/45.23842
2.4184/47.65682

df.mag %>% 
  select(meaningful.unique, mag_admx_group) %>% 
  arrange(mag_admx_group) %>% 
  write_csv("output/admixture/Ash_data/mag_two_clusters.csv", na = "")
