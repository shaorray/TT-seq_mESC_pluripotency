# Rui Shao 
# 2021 Aug

# VariantTools algorithm 
# refer to: https://www.bioconductor.org/help/course-materials/2014/BioC2014/Lawrence_Tutorial.pdf
library(GenomicAlignments)
library(Rsamtools)
library(gmapR)
library(VariantTools)

# mm10 genes
gencode.mm10 <-
  importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.vM25.annotation.gtf")
gencode.mm10$gene_id <- gsub("\\..*", "", gencode.mm10$gene_id)
gencode.exon.mm10 <- gencode.mm10[gencode.mm10$type == "exon"]

# top 200 transcribed genes
top_genes <- rowMeans(tx_L_mat) %>% 
  sort(decreasing = T) %>% 
  names() %>%
  intersect.Vector(gencode.mm10$gene_id) %>% 
  head(200) 

# repetitive elements to mask
repeats <- importRanges("/mnt/0E471D453D8EE463/genomeDir/rmsk/mm10_rmsk_TE.gtf")

# input files
bamfiles <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10", pattern = "bam$", full.names = T)

# sample relative sizes
sf <- sapply(bamfiles, function(x) {
  Rsamtools::idxstatsBam(x) %>% 
    as.data.frame() %>%
    dplyr::filter(grepl("chr", .$seqnames)) %>% 
    dplyr::select(mapped)
}) %>% Reduce(cbind, .) %>% SizeFactorCal() %>% 
  `names<-`(., gsub(".*(LRNA_.*|FRNA_.*).Aligned.*", "\\1", bamfiles))

# single nt variant table template
freq_tab <- expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")) %>%
  dplyr::filter(Var1 != Var2) %>%
  dplyr::mutate(Variant = paste(Var1, Var2, sep = "->")) %>%
  dplyr::mutate(Variant = sort(Variant)) %>%
  dplyr::select(Variant)

# Gmap genome
mouseGmapGenome <- GmapGenome(rtracklayer::FastaFile('/mnt/0E471D453D8EE463/genomeDir/mm10/GRCm38.p5.genome.fa'),
                              create = T)
bpp <- BiocParallel::MulticoreParam(4)

# ------------------------------------------------------------------------------- #
# run
# e.g. c("Smtnl2", "Ccdc14", "Camk2g", "Atf7ip")

gene_res_list <- list()

for (gene in top_genes) {
  param <- TallyVariantsParam(genome = mouseGmapGenome, 
                              mask = repeats, 
                              # which = gencode.exon.mm10[gencode.exon.mm10$gene_name == gene],
                              which = gencode.mm10[gencode.mm10$gene_id == gene & gencode.mm10$type == "gene"],
                              indels = TRUE)
  res_list <- list()
  
  for (i in 1:12) {
    sample_name <- gsub(".*RNA_(.*).Aligned.*", "\\1", bamfiles[i])
    
    bam.index = paste0(bamfiles[i], '.bai')
    if(!file.exists(bam.index)) Rsamtools::indexBam(files = bamfiles[i])
    bam.index = paste0(bamfiles[i + 12], '.bai')
    if(!file.exists(bam.index)) Rsamtools::indexBam(files = bamfiles[i + 12])
    
    # LRNA
    tallies_L <- tallyVariants(bamfiles[i + 12], param, BPPARAM = bpp)
    idx <- which(nchar(tallies_L@ref) == 1 & nchar(tallies_L@alt) == 1 & (tallies_L@refDepth / tallies_L@totalDepth > 0.6))
    tab_L <- table(paste(tallies_L@ref, tallies_L@alt, sep = "->")[idx])
    
    count_L <- data.frame(Type = paste(tallies_L@ref, tallies_L@alt, sep = "->"), 
                          altDepth = tallies_L@altDepth,
                          totalDepth = tallies_L@totalDepth) %>% 
      dplyr::slice(idx) %>%
      dplyr::group_by(Type) %>% 
      dplyr::summarise(convertion = sum(altDepth) / sum(totalDepth),
                       sum_altDepth = sum(altDepth), 
                       sum_totalDepth = sum(totalDepth)) %>% 
      as.data.frame() %>%
      `rownames<-`(.$Type) %>%
      `colnames<-`(paste0(colnames(.), "_L"))
    
    # FRNA
    tallies_F <- tallyVariants(bamfiles[i], param, BPPARAM = bpp)
    idx2 <- which(nchar(tallies_F@ref) == 1 & nchar(tallies_F@alt) == 1 & (tallies_F@refDepth / tallies_F@totalDepth > 0.6))
    tab_F <- table(paste(tallies_F@ref, tallies_F@alt, sep = "->")[idx2])
    
    count_F <- data.frame(Type = paste(tallies_F@ref, tallies_F@alt, sep = "->"), 
                          altDepth = tallies_F@altDepth,
                          totalDepth = tallies_F@totalDepth) %>% 
      dplyr::slice(idx2) %>%
      dplyr::group_by(Type) %>% 
      dplyr::summarise(convertion = sum(altDepth) / sum(totalDepth),
                       sum_altDepth = sum(altDepth), 
                       sum_totalDepth = sum(totalDepth)) %>% 
      as.data.frame() %>%
      `rownames<-`(.$Type) %>%
      `colnames<-`(paste0(colnames(.), "_F"))
    
    tmp <- cbind(freq_tab, 
                 freq_tab_L = as.numeric(tab_L[freq_tab$Variant]),
                 freq_tab_F = as.numeric(tab_F[freq_tab$Variant]),
                 enrich_rate = as.numeric(((tab_L / sum(tab_L)) / (tab_F[names(tab_L)] / sum(tab_F)))[freq_tab$Variant]),
                 count_L[freq_tab$Variant, -1], 
                 count_F[freq_tab$Variant, -1]) %>%
      `rownames<-`(NULL) %>%
      list() %>% 
      `names<-`(sample_name)
    
    res_list <- c(res_list, tmp)
  }
  
  gene_res_list <- c(gene_res_list, `names<-`(list(res_list), gene))
  
}

dat_ER <- sapply(names(gene_res_list), function(gene) {
  idx <- lapply(gene_res_list[[gene]], function(tab) sum(tab[, 3], na.rm = T) > 120) %>% Reduce(c,.)
  if (all(idx)) {
    tmp <- data.frame(Sample = rep(names(gene_res_list[[1]]), each = 12),
                      Variant = lapply(gene_res_list[[gene]], function(tab) tab[, 1]) %>% Reduce(c,.),
                      pos_L = lapply(gene_res_list[[gene]],
                                     function(tab) tab[, 2] %>% "/"(sum(.))
                                     ) %>% Reduce(c,.) ,
                      pos_F = lapply(gene_res_list[[gene]],
                                     function(tab) tab[, 3] %>% "/"(sum(.))
                                     ) %>% Reduce(c,.),
                      sum_altDepth_L = lapply(gene_res_list[[gene]],
                                              function(tab) tab[, 6]
                                              ) %>% Reduce(c,.),
                      sum_totalDepth_L = lapply(gene_res_list[[gene]],
                                              function(tab) tab[, 7]
                      ) %>% Reduce(c,.),
                      gene_id = gene)
    tmp$enrich_pos <- tmp$pos_L / tmp$pos_F
    return(tmp)
  } 
}) %>% Reduce("rbind", .) %>% 
  dplyr::filter(!grepl("SL2i|2i_7d", Sample)) %>%
  dplyr::mutate(Sample = factor(Sample, 
                                levels = c("SL_rep1", "SL_rep2", "SL_rep3",
                                           "2i_2d_rep1", "2i_2d_rep2",
                                           "mTORi_1d_rep1", "mTORi_1d_rep2",
                                           "mTORi_2d_rep1", "mTORi_2d_rep2"))) %>%
  dplyr::mutate(State = factor(gsub("_rep.", "", Sample), levels = c("SL", "2i_2d", "mTORi_1d", "mTORi_2d")))

g1 <- ggplot(dat_ER, aes(x = Variant, y = pos_L, fill = State)) +
  geom_boxplot(outlier.size = 0, notch = T) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  ylim(c(0, 0.3)) +
  xlab("") + ylab("\nRel. variant position (LRNA)") +
  theme_setting + 
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), 
        legend.position = "none")

g2 <- ggplot(dat_ER, aes(x = Variant, y = pos_F, fill = State)) +
  geom_boxplot(outlier.size = 0, notch = T) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  ylim(c(0, 0.3)) +
  xlab("") + ylab("\nRel. variant position (FRNA)") +
  theme_setting + 
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), 
        legend.position = "none")

g3 <- ggplot(dat_ER, aes(x = Variant, y = enrich_pos, fill = State)) +
  geom_hline(yintercept = 1.1, color = add.alpha("red", 0.5), size = 2) +
  geom_hline(yintercept = 1, color = add.alpha("grey50", 0.5), size = 1) +
  geom_boxplot(outlier.size = 0, notch = T) +
  scale_y_continuous(limits = c(0.5, 2), breaks = c(0.5, 1, 1.1, 1.5, 2), labels = c(0.5, 1, 1.1, 1.5, 2)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  xlab("") + ylab("\n Position enrichment (LRNA / FRNA)") + 
  theme_setting + 
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), 
        legend.position = "none")  

cmps <- list( c("SL", "2i_2d"), c("SL", "mTORi_1d"), c("SL", "mTORi_2d") )
g4 <- ggplot(dat_ER[dat_ER$Variant == "T->C", ], aes(x = State, y = enrich_pos, fill = State)) +
  geom_hline(yintercept = 1.1, color = add.alpha("red", 0.5), size = 2) +
  geom_boxplot(outlier.size = 0, notch = T) +
  ggpubr::stat_compare_means(comparisons = cmps, 
                             aes(label = ..p.signif..),
                             method = "t.test", label.y = c(1.7, 1.8, 1.9)) +
  scale_y_continuous(limits = c(0.5, 2), breaks = c(0.5, 1, 1.1, 1.5, 2), labels = c(0.5, 1, 1.1, 1.5, 2)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  xlab("") + ylab("\nT->C enrichment (LRNA / FRNA)") +
  theme_setting + 
  theme(axis.text.x =  element_text(angle = 45, hjust = 1))

g5 <- data.frame(sf = sf[grep("LRNA(_2i_2|_SL_|_mTORi)", names(sf))],
           sample = gsub("LRNA_", "", names(sf[grep("LRNA(_2i_2|_SL_|_mTORi)", names(sf))]))
           ) %>%
  dplyr::mutate(State = factor(gsub("_rep.", "", sample),
                               levels = c("SL", "2i_2d", "mTORi_1d", "mTORi_2d")),
                sample = factor(sample, levels = sample[c(7:9, 1:6)])) %>%
  ggplot(aes(x = sample, y = sf, fill = State)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  xlab("") + ylab("\nRel. LRNA library sizes") +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

g6 <- aggregate(dat_ER[, c("sum_altDepth_L", "sum_totalDepth_L")], 
          by = list(dat_ER$Variant, dat_ER$Sample), "sum") %>%
  as.data.frame() %>% 
  dplyr::mutate(State = factor(gsub("_rep.", "", Group.2), 
                               levels = c("SL", "2i_2d", "mTORi_1d", "mTORi_2d"))
                ) %>%
ggplot(aes(x = Group.1, y = sum_altDepth_L, fill = State)) +
  geom_bar(stat = "identity", width = 0.6,
           position = position_dodge(width = 0.7)) +
  xlab("") + ylab("LRNA variant total depth") +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

ggsave(plot = cowplot::plot_grid(g5, g6, g1, g2, g3, g4, nrow = 3),
       filename = paste0("FigS2_4sU_labeling_comparison.png"), 
       path = "../figS2/figs/",
       device = "png", width = 10, height = 12)
