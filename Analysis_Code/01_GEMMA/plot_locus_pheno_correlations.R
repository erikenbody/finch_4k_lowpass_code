library(noia)
library(tidyverse)

#setup data
df.all.genos <- read.csv("output/combined_genotypes/Darwins_finches_combined_genotype_calls.csv", na.strings = "")

#get tidy names
df_tidy <- read.csv("output/GEMMA/effect_sizes/ld_prune_locus_names_edit.csv")

#update sex with field sex if genetic is missing
df.all.genos <- df.all.genos %>% 
  mutate(genetic.sex = as.factor(gsub("unassigned", NA, genetic.sex)),
         Sex2 = ifelse(Sex == 1, "Male", 
                       ifelse(Sex ==2, "Female", NA))) %>% 
  mutate(munch.sex = ifelse(is.na(genetic.sex) & !is.na(Sex2), as.character(Sex2), as.character(genetic.sex)))


df.dap.genos <- df.all.genos %>% filter(Island == "Daphne")

df.dap.genos <- df.dap.genos 
df.dap.genos[df.dap.genos=="AB"]<-"SL"
df.dap.genos[df.dap.genos=="AA"]<-"SS"
df.dap.genos[df.dap.genos=="BB"]<-"LL"

# get loci of interest ----------------------------------------------------

asm28loci <- grep("asm28loci", names(df.all.genos), value = T)

asm28loci.minus.linked <- grep(paste(c("_2$","_4","_5","_6"), collapse="|"), asm28loci, invert=T, value=T)
asm28loci.minus.gwas <- grep(paste(c("_1$","_3","_7","_27"), collapse="|"), asm28loci.minus.linked, invert=T, value=T)
asm28loci.minus.gwas.with.linked <- grep(paste(c("_1$","_3","_7","_27"), collapse="|"), asm28loci, invert=T, value=T)

all_loci <- c(gwas_loci, asm28loci.minus.gwas)

#some loci highly collinear in scandens
scandens_loci <- all_loci[-c(7,8)]

ld_filt_loc <- c(gwas_loci, "asm28loci_chr2_10", "asm28loci_chr2_12","asm28loci_chr2_17",
                 "asm28loci_chr2_19","asm28loci_chr3_22", "asm28loci_chr3_24", "asm28loci_chr5_25",
                 "asm28loci_chr7_26", "asm28loci_chr25_28")
#gwas_loci <- c("gwas_genotype_chr1_1", "gwas_genotype_chr1A_19", "ALX1.simple",
#               "HMGA2.simple", "gwas_genotype_chr2_21", "gwas_genotype_chr9_23")
#V2
gwas_loci <- c("gwas_genotype_chr1_2", "gwas_genotype_chr1A_17", "ALX1.simple",
               "HMGA2.simple", "gwas_genotype_chr2_18", 
               "gwas_genotype_chr9_20")

# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
# plot all pairwise epistatic effects
##very helpful little function from:
##https://stackoverflow.com/questions/3735286/create-a-matrix-of-scatterplots-pairs-equivalent-in-ggplot2
gatherpairs <- function(data, ..., 
                        xkey = '.xkey', xvalue = '.xvalue',
                        ykey = '.ykey', yvalue = '.yvalue',
                        na.rm = FALSE, convert = FALSE, factor_key = FALSE) {
  vars <- quos(...)
  xkey <- enquo(xkey)
  xvalue <- enquo(xvalue)
  ykey <- enquo(ykey)
  yvalue <- enquo(yvalue)
  
  data %>% {
    cbind(gather(., key = !!xkey, value = !!xvalue, !!!vars,
                 na.rm = na.rm, convert = convert, factor_key = factor_key),
          dplyr::select(., !!!vars)) 
  } %>% gather(., key = !!ykey, value = !!yvalue, !!!vars,
               na.rm = na.rm, convert = convert, factor_key = factor_key)
}

plot_grid_epi <- function(cluster, target_loci, pheno, outname, x, y){

  df.finch.sub <- df.dap.genos %>% 
    #filter(grm_cluster == cluster & !is.na(bill.depth) & !is.na(bill.length) & !is.na(bill.width))
    filter(grm_cluster == cluster) 
  
  df_keys <- df.finch.sub %>% 
    gatherpairs(all_of(target_loci)) %>% 
    group_by(.xkey, .xvalue, .ykey, .yvalue) %>% 
    summarise(phenotype = mean(.data[[pheno]], na.rm = T),
              n = n(), .groups = "keep") %>%
    filter(n > 5) %>% 
    drop_na() 
  
  df_keys2 <- left_join(df_keys, df_tidy[,c("locus", "tidy_locus")], by = c(".xkey" = "locus")) %>% 
    mutate(.xkey = tidy_locus) %>% select(-tidy_locus) %>% 
    left_join(df_tidy[,c("locus", "tidy_locus")], by = c(".ykey" = "locus")) %>% 
    mutate(.ykey = tidy_locus) %>% select(-tidy_locus)
  
  df_keys2$.xvalue <- factor(df_keys2$.xvalue, levels = c("SS", "SL","LL"))
  df_keys2$.yvalue <- factor(df_keys2$.yvalue, levels = c("SS", "SL","LL"))
  
  p1 <- df_keys2 %>% {
      ggplot(., aes(x = .xvalue, y = .yvalue, fill = phenotype)) +
        geom_tile() + 
        geom_text(aes(label = n)) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        facet_wrap(.xkey ~ .ykey, ncol = length(unique(.$.ykey)), scales = 'free', labeller = label_both) +
        theme_bw() +
        ggtitle(paste0("Phenotype = ", pheno))
  }
  
  pdf(paste0("output/pairwise_geno_pheno_correlations/all_pairwise_",outname, "_", pheno,".pdf"), width = x, height = y)
  print(p1)
  dev.off()
}

plot_grid_epi("cluster1",gwas_loci, "bill.PC1.sp", "fortis_gwas_bill.PC1.sp", 18, 16)
plot_grid_epi("cluster1",gwas_loci, "bill.PC2.sp", "fortis_gwas_bill.PC2.sp", 18, 16)

plot_grid_epi("cluster2",gwas_loci, "bill.PC1.sp", "scandens_gwas_bill.PC2.sp", 18, 16)
plot_grid_epi("cluster2",gwas_loci, "bill.PC2.sp", "scandens_gwas_bill.PC2.sp", 18, 16)

plot_grid_epi("cluster1",gwas_loci, "weight", "fortis_gwas_weight", 18, 16)

#plot_grid_epi("fortis", ld_filt_loc, "bill.PC1", "fortis_ld_filt_loc_bill.PC1", 18, 16)
#plot_grid_epi("fortis",all_loci, "bill.PC1", "fortis_all_loci_bill.PC1", 60, 55)
#plot_grid_epi("fortis",all_loci, "bill.PC2", "fortis_all_loci_bill.PC2", 60, 55)


# 3d plot -----------------------------------------------------------------
#not working
species <- "fortis"
pheno <- "bill.PC1"
df.finch.sub <- df.dap.genos.num %>% filter(Species == species & !is.na(bill.depth) & !is.na(bill.length) & !is.na(bill.width))

#only do for the focal subset
beak.pca <- prcomp(df.finch.sub[,c("bill.depth", "bill.width", "bill.length")], center = T, scale = T)
summary(beak.pca)

df.finch.sub <- cbind(df.finch.sub,beak.pca$x[,1:3])
#replace existing values (which was built from 4 species)
df.finch.sub$PC1 <- NULL
df.finch.sub$PC2 <- NULL
df.finch.sub$bill.PC1 <- df.finch.sub$PC1
df.finch.sub$bill.PC2 <- df.finch.sub$PC2
df.finch.sub$bill.PC3 <- df.finch.sub$PC3

names(df.finch.sub) <- gsub("gwas_genotype_", "GW.", names(df.finch.sub))

names(df.finch.sub) <- gsub("asm28loci_", "AS.", names(df.finch.sub))
target_loci<-ld_filt_loc
target_loci <- gsub("gwas_genotype_", "GW.", target_loci)
target_loci <- gsub("asm28loci_", "AS.",target_loci)

df.finch.sub2<-df.finch.sub %>% 
  select(bill.PC1, GW.chr1A_17, GW.chr2_18)

df.mat <- df.finch.sub2 %>% 
  gatherpairs(all_of(c("GW.chr1A_17", "GW.chr2_18"))) %>% 
  group_by(.xkey, .xvalue, .ykey, .yvalue) %>% 
  summarise(phenotype = mean(.data[[pheno]], na.rm = T),
            n = n(), .groups = "keep") %>%
  filter(n > 5) %>% 
  drop_na() 



library(latticeExtra)
#cloud(.xvalue~phenotype +.yvalue, df.mat, panel.3d.cloud=panel.3dbars, col.facet='grey', 
#      xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
#      par.settings = list(axis.line = list(col = "transparent")))


df.mat %>% {
    ggplot(., aes(x = .xvalue, y = .yvalue, fill = phenotype)) +
      geom_tile() + 
      #geom_text(aes(label = round(phenotype,1))) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      facet_wrap(.xkey ~ .ykey, ncol = length(unique(.$.ykey)), scales = 'free', labeller = label_both) +
      theme_bw() +
      ggtitle(paste0("Phenotype = ", pheno))
  }

# linkage among sites -----------------------------------------------------

df.finch.sub <- df.dap.genos.num %>% filter(Species == "fortis")


#for cramers v
library(DescTools)
#for ggplot correlation
library(ggcorrplot)
library(gtools) # combination
library(data.table) # data mgmt

all_loci_linked <- c(gwas_loci, asm28loci.minus.gwas.with.linked)

tmp.df <- df.finch.sub %>% select(all_of(all_loci_linked))
#CramerV(table(tmp.df))


#need to create G prefix names for all loci
name_fixing <- tibble(locus = names(tmp.df)) # %>% write_csv("output/GEMMA/effect_sizes/tmpnames.csv")

nf2 <- full_join(name_fixing, df_tidy, by = "locus") 
nf2 <- nf2 %>% 
  separate(locus, into = c(NA, NA, "suffix"), remove = F) %>% 
  mutate(suffix = as.numeric(ifelse(grepl("chr", suffix), NA, suffix))) %>% 
  mutate(tidy_locus = ifelse(is.na(tidy_locus) & suffix <10, paste0("G0", suffix), 
                             ifelse(is.na(tidy_locus) & suffix > 10, paste0("G", suffix), tidy_locus))) 

#replace the names of tmp.df with the tidy names 
names(tmp.df) <- nf2$tidy_locus[order(match(nf2$locus, names(tmp.df)))]

pdf("output/GEMMA/genotype_plots/correl_genotype_clusters.pdf", width = 30, height = 30)
model.matrix(~0+., data=tmp.df) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)
dev.off() 

ggsave("output/GEMMA/genotype_plots/correl_genotype_clusters.png", width = 6, height = 6)

#cramers V from: https://github.com/Julien-Yacine/cramer_v_heatmap
cv.test = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
              (length(x)[1] * (min(length(unique(x))[1],length(unique(y))[1]) - 1)))
  return(as.numeric(CV))
}

# Apply the function to the combination of categorical variable
v_cramer_all <- function(cat_var, df){
  cat_var_grid <- data.table(combinations(n = length(cat_var), r = 2, v = cat_var, repeats.allowed = T))
  
  do.call(rbind,
          apply(cat_var_grid, 1, function(x){
            tmp <- as.character(x)
            vec1 <- base::unlist(df[,tmp[1]])
            vec2 <- base::unlist(df[,tmp[2]])
            
            data.table(
              variable_x = tmp[1],
              variable_y = tmp[2],
              chi2 = chisq.test(x = vec1, vec2, correct=FALSE)$p.value,
              v_cramer = cv.test(x = vec1, y = vec2)
            )
          }))
  
}

#all_loci_linked.order <- all_loci_linked[c(1,7,4,8:10,3,11:12,2,13:22,5,23:29,6,30)]
#order of tidy names
all_loci_linked.order <- c("G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09","G29", "G10", "G11", "G12", "G13",
"G14", "G15", "G16", "G17", "G18", "G19","G30", "G20", "G21", "G22", "G23", "G24", "G25", "G26",
"G27", "G28")

df.order <- tibble(locus = all_loci_linked.order) %>% 
  mutate(order = 1:n())

results <- v_cramer_all(cat_var = colnames(tmp.df), df = tmp.df)

results$variable_x <- factor(results$variable_x, levels = all_loci_linked.order)
results$variable_y <- factor(results$variable_y, levels = all_loci_linked.order)

results.ax2 <- results %>% 
  dplyr::rename(variable_x = variable_y, variable_y = variable_x)

results2 <- rbind(results, results.ax2)

p.cramers <- ggplot(results2, aes(variable_x, variable_y)) +
  geom_tile(aes(fill = v_cramer), colour = "black") +
  geom_text(aes(label = round(v_cramer, 2)), size = 6) +
  #theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 22),
        axis.text.y=element_text(size = 22)) 
p.cramers

ggsave("output/GEMMA/genotype_plots/cramersV_genotype_clusters_V2.png", width = 20, height = 20)

