# Rui Shao 2020 Dec

# -----------------------------------------------------------------------------------------------------
# This part includes:
#     1. correlations of chromatin features with gene estimated speeds
#     2. gene and ncRNA elongation speed and time
#     3. pausing time estimation and interpretation
# -----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

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
                       Type = c(rep("RNA/DNA", 4), rep("Accessibility", 2), rep("RNA/DNA", 2), 
                                rep("General TF", 5), rep("Pluripotent TF", 6),
                                rep("Enhancer TF", 6), rep("Remodeler", 3), 
                                rep("Heterochromatin", 4), rep("Histone Variant", 3),
                                rep("Histone Acetylation", 6), rep("Histone Methylation", 6),
                                rep("Histone Ubiquitylation", 2))
                       )
data_cor$Type <- factor(data_cor$Type, levels = unique(data_cor$Type))

# plot
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
  xlab("") + ylab("Est. speed correlation")

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
  return(log10(TT_mat + 1) - log10(pol2_mat[, 1] + 1))
}

mRNA.speed <- estimate_elongation_speed(gene.gr)
ncRNA.speed <- estimate_elongation_speed(TU.nc.mm9)
pausing.speed <- estimate_elongation_speed(pause_sites.gr)

# scaling with measured speed
elongation_speed_table <- read.table("../fig4/data/elongation_speed_table.txt")
elongation_speed_table <- elongation_speed_table[#elongation_speed_table$speed > 0.5 &
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
  scale_y_continuous(name = expression('log'[10]*" Est. speed (kb/min)"), 
                     breaks = c((-2):2), labels = 10^c((-2):2)) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20) -1]) + 
  xlab("") + 
  theme_setting


g2 <- ggplot(data_speed_time, aes(x = Type, y = Time, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(name = expression('log'[10]*" Est. elongation time (min)"), 
                     breaks = c((-2):3), labels = 10^c((-2):3), limits = c(-2, 3.5)) +
  coord_cartesian(ylim = c(-2.3, 3.9)) +
  scale_fill_manual(values = colors_20[c(13, 2, 20)]) + 
  xlab("") + 
  theme_setting

ggsave(plot = cowplot::plot_grid(g1, g2, nrow = 2, align = "v"),
       filename = "Fig4_est_elongation_speed_time_TU_types.png", path = "../fig4/figs",
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
    theme(legend.position = "none")
}

# mRNA speed vs uaRNA speed | 0.1676148
g2.1 <- data.frame(x = mRNA.speed[names(coding.promoters)[queryHits(promoter_mtch)], 1], 
           y = ncRNA.speed[subjectHits(promoter_mtch), 1]) %>%
  plot_scatter(.xlab = "Est. mRNA speed (kb/min)\n",
               .ylab = "Est. uaRNA speed (kb/min)",
               r = "0.17")

# mRNA time vs uaRNA time | 0.2087447
g2.2 <- data.frame(x = mRNA.time[names(coding.promoters)[queryHits(promoter_mtch)], 1], 
           y = ncRNA.time[subjectHits(promoter_mtch), 1]) %>%
  plot_scatter(.xlab = "Est. mRNA elongation (min)\n",
               .ylab = "Est. uaRNA elongation (min)", ylim = c(-1.5, 2.5),
               r = "0.21")

# pausing speed vs uaRNA time | 0.2471008
g2.3 <- data.frame(x = pausing.speed[names(coding.promoters)[queryHits(promoter_mtch)], 1], 
           y = ncRNA.speed[subjectHits(promoter_mtch), 1]) %>%
  plot_scatter(.xlab = "Est. pausing speed (kb/min)\n",
               .ylab = "Est. uaRNA speed (kb/min)", ylim = c(-2, 1.7),
               r = "0.25")

# pausing time vs mRNA time | 0.4840085
rm.idx <- which(pausing.time[,1] %in% names(head(sort(table(pausing.time[,1]), T))))
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
       filename = "Fig4_est_elongation_speed_time_correlation.png", path = "../fig4/figs",
       device = "png", width = 10, height = 6.5)

