library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(GenomicRanges)
##this is how I originally set up groupings to run in PIXY
df.hc <- read_xlsx("database/Finch_high_coverage_sample_list.xlsx")

sp.l <- c("fortis", "fuliginosa", "magnirostris","scandens")
df.hc %>% filter(Species_inferred %in% sp.l) %>%
  filter(Island_inferred!="Daphne") %>%
  select(meaningful.unique, Species_inferred) %>% arrange(Species_inferred) %>%
  write_tsv("data/snmf/four_species_for_AIM.txt", col_names = F)

df.hc %>% filter(Species_inferred %in% sp.l) %>%
  filter(Island_inferred!="Daphne") %>%
  select(meaningful.unique, Species_inferred) %>% arrange(Species_inferred) %>%
  select(meaningful.unique) %>%
  write_tsv("data/snmf/four_species_for_AIM_samples.txt", col_names = F)

#select samples to run in snmf
#V2.mdb <- read_excel("/Users/erikenbody/Dropbox/Darwin\'s\ Finch\ shared/Finch_Database_V2.xlsx", guess_max = 1048576)
#df.final.proj <- V2.mdb %>% filter(!is.na(Lowcov_project) | !is.na(High_coverage_project))
#dap.samps <- df.final.proj %>%
#  filter(Island == "Daphne" & Species!="bigbirds") %>%
#  #select(meaningful.unique, Species,Island, Genotype, nest.Breeding.pop) %>%
#  #select(meaningful.unique) %>%
#  distinct(meaningful.unique, .keep_all = T)
#dap.samps %>% select(meaningful.unique) %>% write_tsv("data/snmf/Daphne_four_species.txt", col_names = F)
#dap.samps <- df.final.proj %>%
#  select(meaningful.unique, Species, Island, TubeNum, Genotype, First.year.min, Sex, bill.length, bill.width, bill.depth, weight, sampleID = meaningful.unique)
#dap.samps %>% write_tsv("data/snmf/Daphne_four_species_metadata.txt", col_names = T)

# -------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
pixy.path <- args[1]
outname <- args[2]

dir.create(outname, showWarnings = F)

#pixy.path <- "data/pixy/output_PERSITE/"
pixy.files <- list.files(pixy.path, full.names = T)
#pixy.files <- grep("chr26", pixy.files, value = T, invert = T)

cat_pixy_fst <- do.call("rbind",lapply(pixy.files,
                                       FUN=function(files){
                                         x <- fread(files, skip = 1)
                                         names(x) <- c("pop1",	"pop2",	"CHROM", "BIN_START",	"BIN_END",	"WEIGHTED_FST", "SNPs")
                                         x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                         x
                                       }))

#df.pixy <- read.table(pixy.path, header = T)
#names(df.pixy) <- c("pop1",	"pop2",	"CHROM", "BIN_START",	"BIN_END",	"WEIGHTED_FST", "SNPs")

cat_pixy_fst <- cat_pixy_fst %>% mutate(contrast = paste0("fst_", pop1, "_", pop2))

###notable change - reset negative values to 0
cat_pixy_fst <- cat_pixy_fst %>% mutate(WEIGHTED_FST = ifelse(WEIGHTED_FST < 0 , 0 , WEIGHTED_FST))

cat_pixy_fst.long <- cat_pixy_fst %>%
  select(-pop1, -pop2, -SNPs)  %>%
  pivot_wider(
    names_from = contrast,
    names_sep = "_",
    values_from = WEIGHTED_FST
  )

#four population PBS based on https://doi.org/10.1038/nature13812

pbs_pixy.long <- cat_pixy_fst.long %>%
  mutate(PBS_magnirostris = (-log10(1-fst_fuliginosa_magnirostris) + -log10(1-fst_fortis_magnirostris) + -log10(1-fst_magnirostris_scandens) -
           -log10(1-fst_fortis_fuliginosa) - -log10(1-fst_fuliginosa_scandens)) / 3,
         PBS_fortis = (-log10(1-fst_fortis_fuliginosa) + -log10(1-fst_fortis_scandens) + -log10(1-fst_fortis_magnirostris) -
                         -log10(1-fst_magnirostris_scandens) - -log10(1-fst_fuliginosa_magnirostris)) / 3,
         PBS_fuliginosa = (-log10(1-fst_fortis_fuliginosa) + -log10(1-fst_fuliginosa_magnirostris) + -log10(1-fst_fuliginosa_scandens) -
                         -log10(1-fst_fortis_scandens) - -log10(1-fst_fortis_magnirostris)) / 3,
         PBS_scandens = (-log10(1-fst_fuliginosa_scandens) + -log10(1-fst_fortis_scandens) + -log10(1-fst_magnirostris_scandens) -
                             -log10(1-fst_fortis_fuliginosa) - -log10(1-fst_fuliginosa_magnirostris)) / 3) %>%
  mutate(idx = paste(CHROM, BIN_START, sep = "_"))
#write.csv(pbs_pixy.long, "output/pixy/pbs_four_taxa.csv")
#reset negative PBS to 0
#pbs_pixy.long <- pbs_pixy.long %>%
#  mutate(PBS_magnirostris = ifelse(PBS_magnirostris < 0 , 0 , PBS_magnirostris),
#         PBS_fortis = ifelse(PBS_fortis < 0 , 0 , PBS_fortis),
#         PBS_fuliginosa = ifelse(PBS_fuliginosa < 0 , 0 , PBS_fuliginosa),
#         PBS_scandens = ifelse(PBS_scandens < 0 , 0 , PBS_scandens))


#identify outliers exceeding 95% confidence interval

v.magnirostris <- quantile(pbs_pixy.long$PBS_magnirostris, 0.95, na.rm = T)[[1]]
v.fortis <- quantile(pbs_pixy.long$PBS_fortis, 0.95, na.rm = T)[[1]]
v.fuliginosa <- quantile(pbs_pixy.long$PBS_fuliginosa, 0.95, na.rm = T)[[1]]
v.scandens <- quantile(pbs_pixy.long$PBS_scandens, 0.95, na.rm = T)[[1]]

#write outliers
magnirostris.sites <- pbs_pixy.long %>% filter(PBS_magnirostris > v.magnirostris) %>% select(CHROM, BIN_START, idx)
fortis.sites <- pbs_pixy.long %>% filter(PBS_fortis > v.fortis) %>% select(CHROM, BIN_START, idx)
fuliginosa.sites <- pbs_pixy.long %>% filter(PBS_fuliginosa > v.fuliginosa) %>% select(CHROM, BIN_START, idx)
scandens.sites <- pbs_pixy.long %>% filter(PBS_scandens > v.scandens) %>% select(CHROM, BIN_START, idx)

aim.sites <- rbind(magnirostris.sites, fortis.sites, fuliginosa.sites, scandens.sites) %>%
  distinct(.keep_all = T) #this only keeps one record if there are sites that occue as outliers in multiple contrasts

aim.sites.filt <- rbind(magnirostris.sites, fortis.sites, fuliginosa.sites, scandens.sites) %>%
  group_by(idx) %>% filter(n() == 1) #this drops any site occurring multiple  times

aim.sites %>% ungroup() %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_aims.txt"), col_names = F)

aim.sites.filt %>% ungroup() %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_aims_filt.txt"), col_names = F)


# only aims for subset of species -----------------------------------------

nonfulig <- rbind(magnirostris.sites, fortis.sites, scandens.sites)
nonfulignonscan <- rbind(magnirostris.sites, fortis.sites)

fuliginosa.sites %>%
  filter(!idx %in% nonfulig$idx) %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_fuliginosa_aims.txt"), col_names = F)

rbind(fuliginosa.sites, scandens.sites) %>%
  filter(!idx %in% nonfulignonscan$idx) %>%
  group_by(idx) %>% filter(n() == 1) %>%
  ungroup() %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_scandens_fuliginosa_aims.txt"), col_names = F)

df.sites <- aim.sites.filt

pbs_pixy.long %>%
  filter(idx %in% df.sites$idx) %>%
  select(CHROM, BIN_START, contains('PBS'), POS = BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_aims_with_PBS.txt"), col_names = T)

pbs_pixy.long %>%
  select(CHROM, BIN_START, contains('PBS'), POS = BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_allsites_PBS.txt"), col_names = T)


# exclude selection regions -----------------------------------------------
#retrieve 28 genomic regions identified in finch assembly paper
#important to keep in mind that scandens was not used in that paper to identify regions under selection

df.bed <- read.table("autosomes_for_ful_mag_species_lmm_PEAKS.bed")
names(df.bed) <- c("chr", "start", "stop", "peak")
#extend regions to +/- 50kb
df.bed$start.ext <- df.bed$start - 50000
df.bed$stop.ext <- df.bed$stop + 50000

range.bed <- GRanges(df.bed$chr, IRanges(as.numeric(df.bed$start.ext), as.numeric(df.bed$stop.ext)))
range.pbs <- GRanges(aim.sites.filt$CHROM, IRanges(aim.sites.filt$BIN_START, aim.sites.filt$BIN_START))
comp_overlap <- suppressWarnings(findOverlaps(range.bed, range.pbs))
#comp.annotated <- cbind(range.bed[queryHits(comp_overlap), ], range.pbs[subjectHits(comp_overlap), c(-1, -2, -3)])

pbs.overlap <- range.pbs[subjectHits(comp_overlap), ] %>% as.data.frame() %>%
  mutate(idx = paste(seqnames, start, sep = "_"))

##quick check to see that what I did worked when testing chr25
##aim.sites.filt %>% filter(idx %in% pbs.overlap$idx) %>% nrow()
##aim.sites.filt %>% filter(CHROM == "chr25" & BIN_START > 2439769 & BIN_START < 2712348) %>% nrow()

#write only sites not overlapping with the 28 peaks. Calling these putatively "neutral"
aim.sites.neut <- aim.sites.filt %>% filter(!idx %in% pbs.overlap$idx)
aim.sites.neut %>% ungroup() %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_aims_neutral.txt"), col_names = F)

#remove any chromosome with a selection peak
aim.sites.neut %>%
  filter(!CHROM %in% df.bed$chr) %>%
  ungroup() %>%
  select(CHROM, BIN_START) %>%
  write_tsv(paste0(outname, "/", outname, "_aims_neutral_STRICT.txt"), col_names = F)

# plot all PBS values per chr ---------------------------------------------

for (chr.idx in unique(pbs_pixy.long$CHROM)){

  pbs_pixy.sub <- pbs_pixy.long %>% filter(CHROM == chr.idx)

  p1 <- pbs_pixy.sub %>% ggplot() + geom_point(aes(x = BIN_START/1000000, y = PBS_magnirostris)) +
    geom_hline( yintercept = v.magnirostris, color = "red") +
    labs(x = paste0(outname, " (Mb)"), title = "G. magnirostris") + theme_bw()
  p2 <- pbs_pixy.sub %>% ggplot() + geom_point(aes(x = BIN_START/1000000, y = PBS_fortis)) +
    geom_hline( yintercept = v.fortis, color = "red") +
    labs(x = paste0(outname, " (Mb)"), title = "G. fortis") + theme_bw()
  p3 <- pbs_pixy.sub %>% ggplot() + geom_point(aes(x = BIN_START/1000000, y = PBS_fuliginosa)) +
    geom_hline( yintercept = v.fuliginosa, color = "red") +
    labs(x = paste0(outname, " (Mb)"), title = "G. fuliginosa") + theme_bw()
  p4 <- pbs_pixy.sub %>% ggplot() + geom_point(aes(x = BIN_START/1000000, y = PBS_scandens)) +
    geom_hline( yintercept = v.scandens, color = "red") +
    labs(x = paste0(outname, " (Mb)"), title = "G. scandens") + theme_bw()


  png(paste0(outname, "/", outname,"_",chr.idx, "_PBS.png"), width = 12, height = 10, units = 'in', res = 300)
  print((p1 + p2) / (p3 + p4))
  dev.off()
}

#cat_pixy_fst %>% ggplot() + geom_point(aes(x = BIN_START, y = WEIGHTED_FST)) + facet_wrap(~contrast)


#cat_pixy_fst %>% filter(CHROM == "chr29") %>%
#  ggplot() + geom_point(aes(x = BIN_START/1000000, y = WEIGHTED_FST)) +
#  facet_wrap(~contrast) +
#  theme_bw() + labs(x = "chr29 (Mb)", y = "Fst")
#ggsave("output/pixy/chr29_region_of_interest.png", width = 12, height = 10)
