# Rui Shao 2020 Dec

# -----------------------------------------------------------------------------------------------------
# This part includes:
#     1. correlations of chromatin features with gene estimated speeds
#     2. gene and ncRNA elongation speed and time
#     3. pausing time estimation and interpretation
# -----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")
source("../fig4/F4_scr_get_features_density.r")
# ----------------------------------------------------------------------

# use txdb gene references
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(gene.gr$gene_id),
                       keytype = "ENTREZID",
                       columns = "GENEID")
gene.gr$gene_id <- res$GENEID[match(gene.gr$gene_id, res$ENTREZID)]
names(gene.gr) <- gene.gr$gene_id 

# keep only active genes, longer than 2kb
gene.gr <- gene.gr[width(gene.gr) > 2000 & width(gene.gr) < 1000000]
gene.gr <- gene.gr[seqnames(gene.gr) %in% paste0("chr", c(1:19, "X", "Y"))]
txRPK <- readRDS("../data/txRPK_SL_2i.RData")
gene.gr <- sort(gene.gr[na.omit(match(rownames(txRPK)[rowMeans(txRPK) > 0.2], gene.gr$gene_id))])
names(gene.gr) <- gene.gr$gene_id # 10674 genes

# ---------------------------------------------------------------------

gene.features <- fetch_TU_feature(gene.gr[rownames(est_speed_mat)], is.fast = F)
gene.features[is.na(gene.features) | gene.features < 0] <- 0
gene.features1 <- gene.features[, -1] %>% log1p() %>% trim_quantile() %>% scale() %>% as.data.frame()
est_speed_mat <- readRDS("../fig4/data/est_speed_mat.RData") # log10 scale

data_cor <- data.frame(Correlation = multi_variance_explained(gene.features1, 
                                                              est_speed_mat[,1], 
                                                              is.cor = T),
                       Feature = colnames(gene.features)[-1],
                       Type = c(rep("RNA/DNA", 4), 
                                rep("Accessibility", 2),
                                rep("RNA/DNA", 2), 
                                rep("General TF", 5),
                                rep("Pluripotent TF", 6),
                                rep("Enhancer TF", 6),
                                rep("Remodeler", 3), 
                                rep("Heterochromatin", 4),
                                rep("Histone Variant", 3),
                                rep("Histone Acetylation", 6),
                                rep("Histone Methylation", 6),
                                rep("Histone Ubiquitylation", 2))
                       )
data_cor$Type <- factor(data_cor$Type, levels = unique(data_cor$Type))

# (Fig 4) plot
library(ggpubr)
g1 <- ggdotchart(data_cor, x = "Feature", y = "Correlation",
           color = "Type",                               
           palette = colors_20, 
           sorting = "descending",                     
           add = "segments",                            
           add.params = list(color = "lightgray", size = 1), 
           group = "Type",                             
           ggtheme = theme_pubr()                   
) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray") + 
  xlab("") + ylab("Est. velocity correlation")

ggsave("Fig4_Est_speed_chromatin_feature_explanation.png", 
       plot = g1,
       path = "../fig4/figs", 
       width = 12, height = 5)

# -------------------------------------------------------------------------
# 
pause_sites.gr <- readRDS("../fig4/data/tss_pause_sites.RData")

TU.DE.mm9 <- readRDS("../fig1/data/TU.DE.mm9.RData")
TU.nc.mm9 <- TU.DE.mm9[TU.DE.mm9$location %in% c("intergenic", "uaRNA")]

estimate_elongation_speed <- function(TU.gr) {
  pol2_mat <- .countBam(bam_files = paste0("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/", 
                                           c("P1_Pol5p_SL_CTR_ALL.mm9.fltd.bam", 
                                             "P1_Pol5p_2i_CTR_ALL.mm9.fltd.bam",
                                             "P1_Pol5p_SL_INH_ALL.mm9.fltd.bam")),
                        intervals = TU.gr,
                        stranded = F, 
                        paired.end = "ignore") 
  pol2s5p_input_size_factor <- readRDS("../fig4/data/pol2s5p_input_size_factor.RData")[c(3,1,4)]
  pol2_mat <- sweep(pol2_mat, 2, pol2s5p_input_size_factor, "/")
  
  L_RNA_sf <- readRDS("../data/LRNA.sizefactor.RC.RData")
  TT_mat <- .countBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                             pattern = "LRNA.*(SL_|_2i_2d|mTORi_).*(rep1|rep2).*bam$", full.names = T),
                      intervals = TU.gr,
                      stranded = T,
                      paired.end = "ignore")
  colnames(TT_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1", colnames(TT_mat))
  TT_mat <- sweep(TT_mat, 2, L_RNA_sf[colnames(TT_mat)], "/")
  TT_mat <- sapply(c("SL", "2i", "mTORi"), function(x) rowMeans(TT_mat[, grep(x, colnames(TT_mat))]))
  return(log10(TT_mat + 1) - log10(pol2_mat + 1))
}

mRNA.speed <- estimate_elongation_speed(gene.gr)
ncRNA.speed <- estimate_elongation_speed(TU.nc.mm9)
pausing.speed <- estimate_elongation_speed(pause_sites.gr)

# scaling with measured speed
elongation_speed_table <- read.table("../fig4/data/elongation_speed_table.txt")
elongation_speed_table <- elongation_speed_table[# elongation_speed_table$speed > 0.5 &
                                                   rownames(elongation_speed_table) %in% names(gene.gr), ]

est_speed_residual <- rowMeans(mRNA.speed[rownames(elongation_speed_table), ]) - log10(elongation_speed_table$speed)

# est speed
mRNA.speed <- mRNA.speed - median(est_speed_residual)
ncRNA.speed <- ncRNA.speed - median(est_speed_residual)
pausing.speed <- pausing.speed - median(est_speed_residual)

# est tx time
mRNA.time <- log10(width(gene.gr) / 1e3) - mRNA.speed
ncRNA.time <- log10(width(TU.nc.mm9) / 1e3) - ncRNA.speed
pausing.time <- log10(width(pause_sites.gr) / 1e3) - pausing.speed

# data to plot
data_speed_time <- data.frame(Type = c(rep("mRNA", nrow(mRNA.speed) * 3),
                                       rep(as.character(TU.nc.mm9$location), 3)))

data_speed_time <- cbind(data_speed_time,
                         rbind(reshape::melt(mRNA.speed)[, 2:3], reshape::melt(ncRNA.speed)[, 2:3]),
                         c(reshape::melt(mRNA.time)[, 3], reshape::melt(ncRNA.time)[, 3])
                         )
colnames(data_speed_time) <- c("Type", "Sample", "Speed", "Time")
data_speed_time$Sample <- factor(data_speed_time$Sample, c("SL", "2i", "mTORi"))
data_speed_time$Type <- factor(data_speed_time$Type, c("mRNA", c("intergenic", "uaRNA")))

g1 <- ggplot(data_speed_time, aes(x = Sample, y = Speed, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = "Est. velocity (kb/min)", 
                     breaks = c((-2):2), labels = 10^c((-2):2)) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20) -1]) + 
  xlab("") + 
  theme_setting


g2 <- ggplot(data_speed_time, aes(x = Type, y = Time, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = "Est. elongation time (min)", 
                     breaks = c((-2):3), labels = 10^c((-2):3), limits = c(-2, 3.5)) +
  coord_cartesian(ylim = c(-2.3, 3.9)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20)]) + 
  xlab("") + 
  theme_setting

ggsave(plot = cowplot::plot_grid(g1, g2, nrow = 2, align = "v"),
       filename = "Fig4_est_elongation_speed_time_TU_types.png",
       path = "../fig4/figs",
       device = "png", width = 4.5, height = 6.5)

# p-values
# SL speed
t.test(data_speed_time$Speed[data_speed_time$Sample == "SL" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "SL" & data_speed_time$Type == "intergenic"]) # 0.0161
t.test(data_speed_time$Speed[data_speed_time$Sample == "SL" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "SL" & data_speed_time$Type == "uaRNA"]) # 0.00000000000000022
# 2i speed
t.test(data_speed_time$Speed[data_speed_time$Sample == "2i" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "2i" & data_speed_time$Type == "intergenic"]) # 0.0161
t.test(data_speed_time$Speed[data_speed_time$Sample == "2i" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "2i" & data_speed_time$Type == "uaRNA"]) # 0.00000000000000022
# mTORi speed
t.test(data_speed_time$Speed[data_speed_time$Sample == "mTORi" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "mTORi" & data_speed_time$Type == "intergenic"]) # 0.7622
t.test(data_speed_time$Speed[data_speed_time$Sample == "mTORi" & data_speed_time$Type == "mRNA"],
       data_speed_time$Speed[data_speed_time$Sample == "mTORi" & data_speed_time$Type == "uaRNA"]) # 0.00000000000000022

t.test(data_speed_time$Time[data_speed_time$Type == "mRNA"],
       data_speed_time$Time[data_speed_time$Type == "uaRNA"])

# pausing time by sample
dat_pausing.time <- reshape::melt(pausing.time)
dat_pausing.time$X2 <- factor(dat_pausing.time$X2, c("SL", "2i", "mTORi"))

ggplot(dat_pausing.time, aes(x = X2, y = value, fill = X2)) +
  geom_violin() +
  geom_boxplot(outlier.color = NA, width = 0.3, size = 0.8) +
  xlab("") +
  scale_y_continuous(name = "Est. pausing time (min)", 
                     breaks = c((-2):3), labels = 10^c((-2):3), limits = c(-2.5, 3)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20)]) + 
  theme_setting +
  theme(legend.position = "none")
ggsave(filename = "Fig4_est_pausing_time_by_sample.png", path = "../fig4/figs",
       device = "png", width = 4, height = 4)

# pair uaRNA to gene
gene_intersect <- intersect(names(gene.gr), names(pause_sites.gr))
coding.promoters <- promoters(gene.gr[gene_intersect])
levels(strand(coding.promoters)) <- c("-", "+", "*")

promoter_mtch <- findOverlaps(coding.promoters, TU.nc.mm9)

pausing_index_mat <- readRDS("../fig4/data/pausing_index_SL_2i.RData")

plot_scatter <- function(dat, .xlab, .ylab, xlim = NULL, ylim = NULL, r) {
  dat <- dat[complete.cases(dat) & abs(rowSums(dat)) < Inf, ]
  ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2, 
             label = paste0("     r = ", r, "\nn = ", nrow(dat))) +
    scale_x_continuous(name = .xlab, breaks = -2:2, labels = 10^(-2:2), limits = xlim) +
    scale_y_continuous(name = .ylab, breaks = -2:2, labels = 10^(-2:2), limits = ylim) +
    scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.9) +
    theme_setting +
    theme(legend.position = "none", 
          axis.title = element_text(size= 12))
}

# mRNA velocity vs uaRNA velocity | 0.1676148
g2.1 <- data.frame(x = mRNA.speed[names(coding.promoters)[queryHits(promoter_mtch)], 1], 
           y = ncRNA.speed[subjectHits(promoter_mtch), 1]) %>%
  plot_scatter(.xlab = "Est. mRNA velocity (kb/min)\n",
               .ylab = "Est. uaRNA velocity (kb/min)",
               r = "0.17")

# pausing velocity vs uaRNA time | 0.2471008
g2.2 <- data.frame(x = pausing.speed[names(coding.promoters)[queryHits(promoter_mtch)], 1], 
                   y = ncRNA.speed[subjectHits(promoter_mtch), 1]) %>%
  plot_scatter(.xlab = "Est. pausing win. velocity (kb/min)\n",
               .ylab = "Est. uaRNA velocity (kb/min)", ylim = c(-2, 1.7),
               r = "0.25")

# mRNA velocity vs pausing velocity | 0.4609146
rm.idx <- which(pausing.time[, 1] %in% names(head(sort(table(pausing.time[, 1]), T))))
g2.3 <- data.frame(x = mRNA.speed[-rm.idx, 1], 
           y = pausing.speed[-rm.idx, 1]) %>%
  plot_scatter(.xlab = "Est. mRNA velocity (kb/min)\n",
               .ylab = "Est. pausing win. velocity (kb/min)", xlim = c(-1.8, 2),
               r = "0.461")

# pausing time vs mRNA time | 0.4840085
g2.4 <- data.frame(x = pausing.time[-rm.idx, 1],  y = mRNA.time[-rm.idx, 1]) %>%
  plot_scatter(.xlab = "Est. pausing time (min)\n",
               .ylab = "Est. mRNA elongation (min)",  ylim = c(-1.5, 3.5),
               r = "0.541")

# pausing time vs pausing index | 0.31112
g2.5 <- data.frame(x = pausing.time[, 1],  y = log10(pausing_index_mat[, 1])) %>%
  plot_scatter(.xlab = "Est. pausing time (min)\n",
               .ylab = "Pausing index",
               r = "0.31")

# elongation time vs pausing index | 0.274937
g2.6 <- data.frame(x = mRNA.time[, 1], y = log10(pausing_index_mat[, 1]) ) %>%
  plot_scatter(.xlab = "Est. mRNA elongation (min)\n", xlim = c(-1, 3),
               .ylab = "Pausing index",
               r = "0.28")

ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5, g2.6, nrow = 2),
       filename = "Fig4_est_elongation_speed_time_correlation.png",
       path = "../fig4/figs",
       device = "png", width = 10, height = 6.5)

# ------------------------------------------------------------------------------------ #
# (Fig EV5) log2FC comparison of LRNA, velocity and pausing
LRNA_log2FC <- mcols(TU.DE.mm10.gr)[!is.na(TU.DE.mm10.gr$gene_id),
                                    c("log2FoldChange_LRNA_sp_2i", "log2FoldChange_LRNA_sp_mTORi")]
rownames(LRNA_log2FC) <- TU.DE.mm10.gr$gene_id[!is.na(TU.DE.mm10.gr$gene_id)]

est_speed_log2FC <- log2(10^(est_speed_mat[, 2:3] - est_speed_mat[, 1]))

pausing_index_log2FC <- log2(pausing_index_mat[, c(1,4)] / pausing_index_mat[, 3])


plot_scatter2 <- function(dat, .xlab, .ylab, xlim = NULL, ylim = NULL, r) {
  dat <- dat[complete.cases(dat) & abs(rowSums(dat)) < Inf, ]
  ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    xlab(.xlab) + ylab(.ylab) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2, 
             label = paste0("     r = ", r, "\nn = ", nrow(dat))) +
    scale_color_viridis_c(option = "B", direction = -1, begin = 0.1, end = 0.9) +
    theme_setting +
    theme(legend.position = "none", 
          axis.title = element_text(size= 12))
}

g1 <- data.frame(x = LRNA_log2FC[rownames(est_speed_log2FC), 1], 
           y = est_speed_log2FC[, 1]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Labeled RNA log2FC", 
                .ylab = "Est. velocity log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("2i 2d vs SL")

g2 <- data.frame(x = LRNA_log2FC[rownames(pausing_index_log2FC), 1], 
           y = pausing_index_log2FC[, 1]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Labeled RNA log2FC", 
                .ylab = "Pausing index log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("2i 2d vs SL")

g3 <- data.frame(x = pausing_index_log2FC[rownames(est_speed_log2FC), 1], 
           y = est_speed_log2FC[, 1]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Pausing index log2FC", 
                .ylab = "Est. velocity log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("2i 2d vs SL")

g4 <- data.frame(x = LRNA_log2FC[rownames(est_speed_log2FC), 2], 
           y = est_speed_log2FC[, 2]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Labeled RNA log2FC", 
                .ylab = "Est. velocity log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("mTORi 1d vs SL")

g5 <- data.frame(x = LRNA_log2FC[rownames(pausing_index_log2FC), 2], 
           y = pausing_index_log2FC[, 2]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Labeled RNA log2FC", 
                .ylab = "Pausing index log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("mTORi 1d vs SL")

g6 <- data.frame(x = pausing_index_log2FC[rownames(est_speed_log2FC), 2], 
           y = est_speed_log2FC[, 2]) %>% 
  dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
  plot_scatter2(.xlab = "Pausing index log2FC", 
                .ylab = "Est. velocity log2FC", 
                r = round(cor(.$x, .$y), 3)) + ggtitle("mTORi 1d vs SL")

ggsave(plot = grid.arrange(g1, g2, g3, g4, g5, g6, nrow = 2),
       filename = "FigS4_scatter_synthesis_pausing_speed_log2FC_correlation.png",
       path = "../figS4/figs",
       device = "png", width = 10, height = 6.5)

# ------------------------------------------------------------------------------------ #
# (Fig EV5) faster elongating genes in 2i

faster_gene_2i <- rownames(est_speed_log2FC)[est_speed_log2FC[, 1] > 1.5 & is.finite(est_speed_log2FC[, 1])]
slower_gene_2i <- rownames(est_speed_log2FC)[est_speed_log2FC[, 1] < (-3) & is.finite(est_speed_log2FC[, 1])]

enrichGeneSets(gene_id = faster_gene_2i, method = "GO",
               ontology = "BP", is.GeneRatio = F, top_n_term = 15,
               title = "2i faster elongation genes (n = 579, log2FC > 1.5)")
ggsave(filename = "FigS4_GO_BP_2i_faster_speed_gene_pval.png",
       path = "../figS4/figs",
       device = "png", width = 9, height = 8)

enrichGeneSets(gene_id = slower_gene_2i, method = "GO", 
               ontology = "BP", is.GeneRatio = F, top_n_term = 15,
               title = "2i slower elongation genes (n = 449, log2FC < -3)")
ggsave(filename = "FigS4_GO_BP_2i_slower_speed_gene_pval.png",
       path = "../figS4/figs",
       device = "png", width = 9, height = 8)

# export estimated velocity
tmp_m <- as.data.frame(10^mRNA.speed)
colnames(tmp_m) <- paste0(colnames(tmp_m), "_kb_min")
tmp_m <- cbind(gene.gr, tmp_m)

TU.nc.mm9_2 <- TU.nc.mm9
mcols(TU.nc.mm9_2) <- data.frame(gene_id = TU.nc.mm9$location)
tmp_nc <- as.data.frame(10^ncRNA.speed)
colnames(tmp_nc) <- paste0(colnames(tmp_nc), "_kb_min")
tmp_nc <- cbind(TU.nc.mm9_2, tmp_nc)

write.table(rbind(tmp_m, tmp_nc), 
            file = "../fig4/data/Est_transcription_velocity_mm9.txt",
            quote = F, row.names = F)

library(openxlsx)
write.xlsx(rbind(tmp_m, tmp_nc), 
           file = "../fig4/data/Est_transcription_velocity_mm9.xlsx")

# ------------------------------------------------------------------------------ #
# (Author Response to Reviewer 3) velocity neighboring effect 
# the same process as "Fig2_TU_coexpression.R"
if (T) {
  TU.intergenic.mm9.gr <- TU.nc.mm9[TU.nc.mm9$location == "intergenic"]
  TU.intergenic.mm9.gr$id <- seq_len(length(TU.intergenic.mm9.gr))
  
  F5_enhancer_mm9 <- importRanges("../../data/anno_ref/FANTOM5/enhancers/mouse_permissive_enhancers_phase_1_and_2_mm9.bed")
  
  TU.intergenic.mm9_as <- promoters(TU.intergenic.mm9.gr, upstream = 500, downstream = 0)
  levels(strand(TU.intergenic.mm9_as)) <- c("-", "+", "*")
  
  mtch_dir <- findOverlaps(TU.intergenic.mm9.gr, TU.intergenic.mm9_as)
  TU.intergenic.mm9.gr$direction <- ifelse(countQueryHits(mtch_dir) > 0,
                                       "Bidirectional", "Unidirectional" )
  TU.intergenic.mm9.gr$enhancer <- ifelse(findOverlaps(F5_enhancer_mm9, TU.intergenic.mm9.gr) %>% countSubjectHits > 0,
                                      "Enhancer", "TX")
  
  
  gene.intergenic.pairs.mm9 <- 
    foreach (i = seqlevels(TU.intergenic.mm9.gr), .combine = rbind) %dopar% {
      # extract intergenic neighbors for each chromosome
      genes.tmp = gene.gr[seqnames(gene.gr) == i]
      genes.tmp = genes.tmp[order(start(genes.tmp))]
      TU.tmp = TU.intergenic.mm9.gr[seqnames(TU.intergenic.mm9.gr) == i]
      TU.strands = as.character(strand(TU.tmp)) == "+"
      
      out_table <- NULL
      for (j in seq_along(genes.tmp)) {
        gene.strand = as.character(strand(genes.tmp[j])) == "+"
        gaps = start(TU.tmp) - start(genes.tmp[j])
        
        same.strand = (TU.strands == gene.strand)
        
        if (j == 1) {
          up.tu =  .ifelse(gene.strand, gaps < 0, start(TU.tmp) < start(genes.tmp[2]) & gaps > 0)
          down.tu = .ifelse(!gene.strand, gaps < 0, start(TU.tmp) < start(genes.tmp[2]) & gaps > 0)
        } else if (j == length(genes.tmp)) {
          up.tu =  .ifelse(gene.strand, 
                           start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                           gaps > 0)
          down.tu = .ifelse(!gene.strand,
                            start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                            gaps > 0)
        } else {
          up.tu =  .ifelse(gene.strand, 
                           gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]), 
                           gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
          down.tu = .ifelse(!gene.strand, 
                            gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]),
                            gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
        }
        up.s = same.strand & up.tu
        up.as = !same.strand & up.tu
        down.s = same.strand & down.tu
        down.as = !same.strand & down.tu
        
        out_table <- rbind(out_table,
                           data.frame("gene_id" = rep(genes.tmp[j]$gene_id, sum(up.s + up.as + down.s + down.as)),
                                      "strand" = rep(strand(genes.tmp[j]), sum(up.s + up.as + down.s + down.as)),
                                      "type" = c(rep("up.sense", sum(up.s)), 
                                                 rep("up.antisense", sum(up.as)),
                                                 rep("down.sense", sum(down.s)),
                                                 rep("down.antisense", sum(down.as)) ),
                                      "gap" = c(gaps[up.s], gaps[up.as], gaps[down.s], gaps[down.as]),
                                      "id" = c(TU.tmp$id[up.s], TU.tmp$id[up.as], TU.tmp$id[down.s], TU.tmp$id[down.as])))
      }
      out_table
    }
  
  # gap of TU TSS to gene boundary
  names(gene.gr) <- gene.gr$gene_id
  gaps <- start(promoters(TU.intergenic.mm9.gr[gene.intergenic.pairs.mm9$id], upstream = 0, downstream = 0)) - 
    cbind(start(gene.gr[gene.intergenic.pairs.mm9$gene_id]),
          end(gene.gr[gene.intergenic.pairs.mm9$gene_id]))
  idx <- ifelse(gene.intergenic.pairs.mm9$strand == "+", 
                ifelse(grepl("up", gene.intergenic.pairs.mm9$type), 1, 2), 
                ifelse(grepl("up", gene.intergenic.pairs.mm9$type), 2, 1))
  
  gene.intergenic.pairs.mm9$gap_boundary <- 
    sapply(seq_along(idx),
           function(x) gaps[x, idx[x]]
    ) * ifelse(gene.intergenic.pairs.mm9$strand == "+", 1, -1) / 1e3
  
  tmp_table <- cbind(mRNA.speed[gene.intergenic.pairs.mm9$gene_id, ], 
                     ncRNA.speed[TU.nc.mm9$location == "intergenic", ][gene.intergenic.pairs.mm9$id, ])
  colnames(tmp_table) <- paste0(c(rep("Gene_", 3), rep("TU_", 3)), colnames(tmp_table))
  gene.intergenic.pairs.mm9 <- cbind(gene.intergenic.pairs.mm9, tmp_table)
  gene.intergenic.pairs.mm9$enhancer <- TU.intergenic.mm9.gr$enhancer[gene.intergenic.pairs.mm9$id]
  gene.intergenic.pairs.mm9$direction <- TU.intergenic.mm9.gr$direction[gene.intergenic.pairs.mm9$id]
  rm(tmp_table)
  
  gap_breaks <- c(min(gene.intergenic.pairs.mm9$gap_boundary), 
                  c(-400, -200, -100, -75, -50, -30, -20, -10, 0, 10, 20, 30, 50, 75, 100, 200, 400),# * 1000,
                  max(gene.intergenic.pairs.mm9$gap_boundary))
  
  gene.intergenic.pairs.mm9$gap_class <- cut(gene.intergenic.pairs.mm9$gap_boundary, 
                                             breaks = gap_breaks, 
                                             labels = 1:18) %>% as.numeric()
  gene.intergenic.pairs.mm9 <- gene.intergenic.pairs.mm9[!is.na(gene.intergenic.pairs.mm9$gap_class), ]

  # TU ~ gene velocity correlation
  mat <- NULL
  for( i in 1:18 ) {
    tmp.pairs <- gene.intergenic.pairs.mm9[as.numeric(gene.intergenic.pairs.mm9$gap_class) == i, ]
    as.idx <- grepl("antisense", tmp.pairs$type)
    
    enh.idx <- grepl("Enhancer", tmp.pairs$enhancer)
    bi.idx <- grepl("Bidirectional", tmp.pairs$direction)
    
    mat <- rbind(mat, c("Pos" = i, 
                        "cor_SL_Sense" = cor(tmp.pairs[!as.idx, "Gene_SL"], tmp.pairs[!as.idx, "TU_SL"]),
                        "cor_2i_Sense" = cor(tmp.pairs[!as.idx, "Gene_2i"], tmp.pairs[!as.idx, "TU_2i"]),
                        "cor_mTORi_Sense" = cor(tmp.pairs[!as.idx, "Gene_mTORi"], tmp.pairs[!as.idx, "TU_mTORi"]),
                        
                        "cor_SL_Antisense" = cor(tmp.pairs[as.idx, "Gene_SL"], tmp.pairs[as.idx, "TU_SL"]),
                        "cor_2i_Antisense" = cor(tmp.pairs[as.idx, "Gene_SL"], tmp.pairs[as.idx, "TU_2i"]),
                        "cor_mTORi_Antisense" = cor(tmp.pairs[as.idx, "Gene_mTORi"], tmp.pairs[as.idx, "TU_mTORi"]),
                        
                        "cor_SL_Enhancer" = cor(tmp.pairs[enh.idx, "Gene_SL"], tmp.pairs[enh.idx, "TU_SL"]),
                        "cor_2i_Enhancer" = cor(tmp.pairs[enh.idx, "Gene_SL"], tmp.pairs[enh.idx, "TU_2i"]),
                        "cor_mTORi_Enhancer" = cor(tmp.pairs[enh.idx, "Gene_mTORi"], tmp.pairs[enh.idx, "TU_mTORi"]),
                        
                        "cor_SL_Non_enh_Bidirectional" = cor(tmp.pairs[!enh.idx&bi.idx, "Gene_SL"], tmp.pairs[!enh.idx&bi.idx, "TU_SL"]),
                        "cor_2i_Non_enh_Bidirectional" = cor(tmp.pairs[!enh.idx&bi.idx, "Gene_SL"], tmp.pairs[!enh.idx&bi.idx, "TU_2i"]),
                        "cor_mTORi_Non_enh_Bidirectional" = cor(tmp.pairs[!enh.idx&bi.idx, "Gene_mTORi"], tmp.pairs[!enh.idx&bi.idx, "TU_mTORi"])
    ))
  }
  mat[10:18, 1] <- (10:18) + 3 # add gene box
  mat <- rbind(mat[1:9, ], matrix(c(10:12, rep(NA, 36)), nrow = 3), mat[10:18, ])
  
  dat_mat <- reshape2::melt(data = as.data.frame(mat), id = "Pos")
  dat_mat$Sample <- factor(gsub(".*(SL|2i|mTORi).*", "\\1", dat_mat$variable), 
                           levels = c("SL", "2i", "mTORi"))
  dat_mat$Direction <- factor(gsub("cor_(.*)_(.*)", "\\2", dat_mat$variable), 
                              levels = c("Sense", "Antisense", "Enhancer", "Bidirectional"))
  dat_mat$line_group <- factor(cumsum(is.na(dat_mat$value)))
  
  ggplot(dat_mat[dat_mat$Direction %in% c("Sense", "Antisense"), ],
         aes(x = Pos, y = value, color = Direction, group = line_group)) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = -Inf, ymax = 0),
              fill = "#0000B2", linetype = 0) +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = 0, ymax = Inf),
              fill = "grey75", linetype = 0) +
    geom_rect(aes(xmin = 9, xmax = 10, ymin = -Inf, ymax = Inf),
              fill = "grey95", linetype = 0) +
    geom_rect(aes(xmin = 12, xmax = 13, ymin = -Inf, ymax = Inf),
              fill = "grey95", linetype = 0) +
    geom_point() + 

    geom_line() +
    
    facet_grid(Sample~.) + 
    scale_x_continuous(name = "\nDistance to gene (kb)", 
                       breaks = c(1, 5, 8, 10, 12, 14, 17, 21), 
                       labels = c(gap_breaks[c(2, 5, 8)], c("5\'", "3\'"), gap_breaks[c(12, 15, 18)]) ) +
    scale_y_continuous(breaks = c(-0.2, 0, 0.2),
                       labels = c(-0.2, 0, 0.2), 
                       limits = c(-0.3, 0.3)) +
    scale_color_manual(name = "Intergenic", values = colors_9[7:6]) +
    ylab("(log10) Est. velocity correlation") +
    theme_setting + 
    theme(strip.text = element_text(size = 12),
          axis.ticks = element_line(), 
          panel.grid.major = element_line(), 
          panel.spacing = unit(1, "lines"),
          legend.position = "top")  
    
  ggsave(filename = "FigS4_est_velocity_intergenic_TU_gene_correlation.png",
         path = "../figS4/figs/", device = "png", width = 5, height = 4)
  
  
  ggplot(dat_mat[dat_mat$Direction %ni% c("Sense", "Antisense"), ],
         aes(x = Pos, y = value, color = Direction, group = line_group)) +
    geom_hline(yintercept = 0, lty = 2, col = "grey50") +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = -Inf, ymax = 0),
              fill = "#0000B2", linetype = 0) +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = 0, ymax = Inf),
              fill = "grey75", linetype = 0) +
    geom_rect(aes(xmin = 9, xmax = 10, ymin = -Inf, ymax = Inf),
              fill = "grey95", linetype = 0) +
    geom_rect(aes(xmin = 12, xmax = 13, ymin = -Inf, ymax = Inf),
              fill = "grey95", linetype = 0) +
    geom_point() + 
    
    geom_line() +
    
    facet_grid(Sample~. ) + 
    scale_x_continuous(name = "\nDistance to gene (kb)", 
                       breaks = c(1, 5, 8, 10, 12, 14, 17, 21), 
                       labels = c(gap_breaks[c(2, 5, 8)], c("5\'", "3\'"), gap_breaks[c(12, 15, 18)]) ) +
    scale_y_continuous(breaks = c(-0.2, 0, 0.2),
                       labels = c(-0.2, 0, 0.2), 
                       limits = c(-0.3, 0.3)) +
    scale_color_manual(name = "Direction", values = colors_9[7:6]) +
    ylab("(log10) Est. velocity correlation") +
    theme_setting + 
    theme(strip.text = element_text(size = 12),
          axis.ticks = element_line(), 
          panel.grid.major = element_line(), 
          panel.spacing = unit(1, "lines"),
          legend.position = "top")  
  
  ggsave(filename = "FigS4_est_velocity_intergenic_TU_gene_enh_bi_correlation.png",
         path = "../figS4/figs/", device = "png", width = 5, height = 4)
  
}




