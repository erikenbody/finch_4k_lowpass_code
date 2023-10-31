library(tidyverse)
library(patchwork)
library(broom)
fort.color <- '#0072B2'
magn.color <- '#E69F00'
scan.color <- '#009E73'
fuli.color <- '#56B4E9'
hybr.color <- '#8da0cb'

df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

df.dap <- df.all.genos %>% filter(Island == "Daphne") %>% 
  mutate(beak_size = bill.length + bill.depth + bill.width,
         beak_shape = bill.length/bill.depth)


df_set <- df.dap %>% 
  filter(Species == "fortis" | Species == "scandens" | 
           Species == "hybrid" & !grepl("f", Genotype)) %>% 
  mutate(G07 = case_when(
    ALX1.simple == "AA" ~ "PP",
    ALX1.simple == "AB" ~ "BP",
    ALX1.simple == "BB" ~ "BB"
    )) 

model <- cor.test(df_set$anc_fuli, df_set$beak_size)

# Get the summary statistics
summary_stats <- glance(model)

# Get R-squared and p-value
r_squared <- summary_stats$estimate ^2
p_value <- summary_stats$p.value


pS <- df_set %>% ggplot() + 
  geom_point(aes(x = bill.length + bill.depth + bill.width, 
                 y = anc_scan, color = G07), alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", aes(x = bill.length + bill.depth + bill.width, y = anc_fuli)) +
  annotate("text", x = 27, y = 0.9, 
           label = paste("r^2 = ", round(r_squared, 4), 
                         "\np-value = ", round(p_value, 2)), hjust = 1.1, vjust = 2) +
  theme_bw() +
  labs(y = "G. scandens ancestry", x = "Bill length / Bill depth") +
  theme(legend.position = "bottom")

df_set2 <- df.dap %>% 
  filter(Species == "fortis" | Species == "fuliginosa" | 
           Species == "hybrid" & !grepl("S", Genotype)) %>% 
  mutate(G03 = case_when(
    HMGA2.simple == "AA" ~ "SS",
    HMGA2.simple == "AB" ~ "SL",
    HMGA2.simple == "BB" ~ "LL"
  ))



model <- cor.test(df_set2$anc_fuli, df_set2$beak_size)

# Get the summary statistics
summary_stats <- glance(model)

# Get R-squared and p-value
r_squared <- summary_stats$estimate ^2
p_value <- summary_stats$p.value

# Add the regression line, R-squared and p-value
pf <- df_set2 %>% 
  ggplot() + 
  geom_point(aes(x = bill.length + bill.depth + bill.width, 
                 y = anc_fuli, color = G03), alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", aes(x = bill.length + bill.depth + bill.width, y = anc_fuli)) +
  annotate("text", x = Inf, y = Inf, 
           label = paste("r^2 = ", round(r_squared, 4), 
                         "\np-value = ", round(p_value, 2)), hjust = 1.1, vjust = 2) +
  theme_bw() +
  labs(y = "G. fuliginosa ancestry", x = "Bill length + Bill depth + Bill width") +
  theme(legend.position = "bottom")

pS + pf
ggsave("output/admixture/ancestry_vs_morpho_g03_g07.png", width =8, height = 6)



