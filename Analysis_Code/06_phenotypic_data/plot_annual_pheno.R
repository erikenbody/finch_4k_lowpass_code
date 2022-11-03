library(tidyverse)
library(patchwork)

fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

head(df.all.genos)

df.pheno.an <- df.all.genos %>% 
  select(meaningful.unique, Species, First.year.min, Last.year, starts_with("bill"), weight)

df.pheno.yr <- df.pheno.an %>% 
  filter(Species!="hybrid") %>% 
  #select(-Species) %>% 
  filter(!is.na(First.year.min)) %>% filter(!is.na(Last.year)) %>% 
  mutate(year = map2(First.year.min, Last.year, `:`)) %>% 
  select(-First.year.min, -Last.year) %>% 
  unnest(cols = c(year))

df.pheno.sum.yr <- df.pheno.yr %>% 
  group_by(year, Species) %>% 
  filter(!is.na(bill.PC1.sp)) %>% 
  summarise_at(.vars = c("bill.PC1.sp", "bill.PC2.sp", "bill.depth", "bill.length", "bill.width", "weight"),
               .funs = list(mean = mean, sd = sd, "se" = ~sd(.x/sqrt(length(.x)))))

write_csv(df.pheno.sum.yr, "output/phenotypic_data/annual_plots/annual_phenotype_summaries.csv")

# plot --------------------------------------------------------------------

df.pheno.sum.yr$Species <- factor(df.pheno.sum.yr$Species, levels = c("fortis", "scandens", "magnirostris"))

p.anPC1 <- df.pheno.sum.yr %>% 
  mutate(lower = bill.PC1.sp_mean - bill.PC1.sp_se, upper = bill.PC1.sp_mean + bill.PC1.sp_se) %>% 
  filter(year > 1982 & Species == "fortis" | 
           year > 1982 & Species == "scandens" | 
           year > 1982 & Species == "magnirostris") %>% 
  ggplot(aes(x = year, y = bill.PC1.sp_mean, color = Species)) +
  geom_point() +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
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
  ylab("Beak PC1") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color))

p.anPC1
ggsave("output/phenotypic_data/annual_plots/beak_pc1_annual.png", width = 8.5, height = 6)
ggsave("output/phenotypic_data/annual_plots/beak_pc1_annual_skinny.png", width = 8.5, height = 3)

p.anPC2 <- df.pheno.sum.yr %>% 
  mutate(lower = bill.PC2.sp_mean - bill.PC2.sp_se, upper = bill.PC2.sp_mean + bill.PC2.sp_se) %>% 
  filter(year > 1982 & Species == "fortis" | 
           year > 1982 & Species == "scandens" | 
           year > 1982 & Species == "magnirostris") %>%   ggplot(aes(x = year, y = bill.PC2.sp_mean, color = Species)) +
  geom_point() +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
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
  ylab("Beak PC2") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color))
p.anPC2
ggsave("output/phenotypic_data/annual_plots/beak_PC2_annual.png", width = 8.5, height = 6)
ggsave("output/phenotypic_data/annual_plots/beak_PC2_annual_skinny.png", width = 8.5, height = 3)


p.anWeight <- df.pheno.sum.yr %>% 
  mutate(lower = weight_mean - weight_se, upper = weight_mean + weight_se) %>% 
  filter(year > 1982 & Species == "fortis" | 
           year > 1982 & Species == "scandens" | 
           year > 1982 & Species == "magnirostris") %>%   ggplot(aes(x = year, y = weight_mean, color = Species)) +
  geom_point() +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
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
  ylab("Weight") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color))
p.anWeight

p.anPC1 <- p.anPC1 + 
  theme(axis.text.x = element_blank())
p.anPC2 <- p.anPC2 + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

p.anPC1 / p.anPC2

ggsave("output/phenotypic_data/annual_plots/beak_PC1_PC2_annual.png", width = 12, height = 5.5)


# just fortis and scan ----------------------------------------------------

p.anPC1.2sp <- df.pheno.sum.yr %>% 
  mutate(lower = bill.PC1.sp_mean - bill.PC1.sp_se, upper = bill.PC1.sp_mean + bill.PC1.sp_se) %>% 
  filter(year > 1982 & Species == "fortis" | 
           year > 1982 & Species == "scandens") %>% 
  ggplot(aes(x = year, y = bill.PC1.sp_mean, color = Species)) +
  geom_point() +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  scale_x_continuous(limits = c(1982.5, 2014.5), breaks = seq(1985, 2010, by = 5)) +
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
  ylab("Beak PC1") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color))

p.anPC2.2sp <- df.pheno.sum.yr %>% 
  mutate(lower = bill.PC2.sp_mean - bill.PC2.sp_se, upper = bill.PC2.sp_mean + bill.PC2.sp_se) %>% 
  filter(year > 1982 & Species == "fortis" | 
           year > 1982 & Species == "scandens") %>%   ggplot(aes(x = year, y = bill.PC2.sp_mean, color = Species)) +
  geom_point() +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  scale_x_continuous(limits = c(1982.5, 2014.5), breaks = seq(1985, 2010, by = 5)) +
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
  ylab("Beak PC2") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color, magn.color))

p.anPC1.2sp <- p.anPC1.2sp + 
  theme(axis.text.x = element_blank())
p.anPC2.2sp <- p.anPC2.2sp + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

p.anPC1.2sp / p.anPC2.2sp

ggsave("output/phenotypic_data/annual_plots/beak_PC1_PC2_annual_2sp.png", width = 8.5, height = 3.5)


# the 3 traits ------------------------------------------------------------



df.pheno.sum.yr %>% 
  mutate(lower = bill.depth_mean - bill.depth_se, upper = bill.depth_mean + bill.depth_se) %>% 
  filter(year > 1982 & Species == "fortis" | year > 1982 & Species == "scandens") %>% 
  ggplot(aes(x = year, y = bill.depth_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak depth") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color))

ggsave("output/phenotypic_data/annual_plots/beak_depth_annual.png", width = 8.5, height = 6)
ggsave("output/phenotypic_data/annual_plots/beak_depth_annual_skinny.png", width = 8.5, height = 3)

df.pheno.sum.yr %>% 
  mutate(lower = bill.length_mean - bill.length_se, upper = bill.length_mean + bill.length_se) %>% 
  filter(year > 1982 & Species == "fortis" | year > 1982 & Species == "scandens") %>% 
  ggplot(aes(x = year, y = bill.length_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak length") +
  facet_grid(~Species) +
  scale_color_manual(values = c(fort.color, scan.color))

ggsave("output/phenotypic_data/annual_plots/beak_length_annual.png", width = 8.5, height = 6)
ggsave("output/phenotypic_data/annual_plots/beak_length_annual_skinny.png", width = 8.5, height = 3)


# species 1 by 1 ----------------------------------------------------------

df.pheno.sum.yr %>% 
  mutate(lower = bill.PC1.sp_mean - bill.PC1.sp_se, upper = bill.PC1.sp_mean + bill.PC1.sp_se) %>% 
  filter(year > 1982 & Species == "fortis") %>% 
  ggplot(aes(x = year, y = bill.PC1.sp_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak Size") +
  #facet_grid(~Species) +
  scale_color_manual(values = c(fort.color))
ggsave("output/phenotypic_data/annual_plots/fortis_beak_pc1_annual_skinny.png", width = 5, height = 3)


df.pheno.sum.yr %>% 
  mutate(lower = bill.PC1.sp_mean - bill.PC1.sp_se, upper = bill.PC1.sp_mean + bill.PC1.sp_se) %>% 
  filter(year > 1982 & Species == "fortis") %>% 
  ggplot(aes(x = year, y = bill.PC1.sp_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak Size") +
  #facet_grid(~Species) +
  scale_color_manual(values = c(fort.color))
ggsave("output/phenotypic_data/annual_plots/just_fortis_pheno.pdf", width = 5, height = 3)

df.pheno.sum.yr %>% 
  mutate(lower = bill.PC2.sp_mean - bill.PC2.sp_se, upper = bill.PC2.sp_mean + bill.PC2.sp_se) %>% 
  filter(year > 1982 & Species == "scandens") %>% 
  ggplot(aes(x = year, y = bill.PC2.sp_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak Shape") +
  #facet_grid(~Species) +
  scale_color_manual(values = c(scan.color))
ggsave("output/phenotypic_data/annual_plots/scandens_beak_pc1_annual_skinny.png", width = 5, height = 3)


df.pheno.sum.yr %>% 
  mutate(lower = bill.PC1.sp_mean - bill.PC1.sp_se, upper = bill.PC1.sp_mean + bill.PC1.sp_se) %>% 
  filter(year > 1982 & Species == "magnirostris") %>% 
  ggplot(aes(x = year, y = bill.PC1.sp_mean, color = Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), width = 0.2) +
  #geom_smooth(size = 2, alpha = 0.6, method = "gam", 
  #            formula = y ~ s(x, bs = "cs", fx = TRUE, k = 22), se = F) +
  scale_x_continuous(limits = c(1982.5, 2012.5), breaks = seq(1985, 2010, by = 5)) +
  #scale_x_continuous(limits = c(1987.5, 2011.5), breaks = seq(1988, 2012, by = 1)) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "italic")) +
  ylab("Beak Size") +
  #facet_grid(~Species) +
  scale_color_manual(values = c(magn.color))
ggsave("output/phenotypic_data/annual_plots/magnirostris_beak_pc1_annual_skinny.png", width = 5, height = 3)


# -------------------------------------------------------------------------


