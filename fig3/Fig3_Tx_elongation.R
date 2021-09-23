#-----------------------------------------------------------------------------------------------------
# This part estimates transcription elongation dynamics with TT-seq and Pol II S5p coverage
# evaluates RNA Pol II pausing in different mouse ES pluripotent states
#
# Rui Shao, June 2020
#-----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

#-----------------------------------------------------------------------------------------------------
# Estimate speeds with TT-seq and Pol II S5p
# Define parameters:
#   v_bar: average elongation speed
#   l: gene length
#   P_0: number of Pol II starting elongation per unit of time
#   P_1: total number of Pol II on gene body engaging on elongation
#   P_bar: number of Pol II per Kb
#   S_bar: RNA synthesis rate (~ TT-seq LRNA RPK or Copy * labeled rate)
#   v_hat: estimated mean speed 
#   m: minutes of elongation

# Equations 
# 1) P_0 ~ S_bar (~ links with a scaling factor)
# 2) P_1 = P_0 * l / v_bar
# 3) local Pol II speed v_i = P0 / d_i, where d_i is local Pol II RPK density
# 4) v_bar = \sum_i^m v_i / m
# 5) P_bar = P_1 / l = P_0 / v_bar


# Pol II elongation speed gene body coverage ------------------------------------------------------------
source("F3_src_Tx_elongation.R")

# load processed data
pol2s5p_sandwich_mat_list <- readRDS("data/pol2s5p_gene_sandwich_mat.RData")
TT_gene_sandwich_mat_list <- readRDS("data/TT_gene_sandwich_mat_norm.RData")
TT_Pol2s5p_sandwich_mat_list <- readRDS("data/TT_Pol2s5p_sandwich_mat.RData") # log1p transformed

# plot average coverage
n_pos <- ncol(pol2s5p_sandwich_mat_list[[1]])

pol2s5p_data <- data.frame(position = rep(seq_len(n_pos), 3),
                           value = c(colMeans(log1p(pol2s5p_sandwich_mat_list[[1]])),
                                     colMeans(log1p(pol2s5p_sandwich_mat_list[[2]])),
                                     colMeans(log1p(pol2s5p_sandwich_mat_list[[3]]))),
                           sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                           levels = c("SL", "2i", "mTORi"))
)

ttseq_data <- data.frame(position = rep(seq_len(n_pos), 3),
                         value = c(colMeans(log1p(TT_gene_sandwich_mat_list[[1]])),
                                   colMeans(log1p(TT_gene_sandwich_mat_list[[2]])),
                                   colMeans(log1p(TT_gene_sandwich_mat_list[[3]]))),
                         sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                         levels = c("SL", "2i", "mTORi"))
)
n_pos <- ncol(pol2s5p_sandwich_mat_list[[1]])
TT_Pol2s5p_data <- data.frame(position = rep(seq_len(n_pos), 3),
                              value = c(colMeans(TT_Pol2s5p_sandwich_mat_list[[1]]),
                                        colMeans(TT_Pol2s5p_sandwich_mat_list[[2]]),
                                        colMeans(TT_Pol2s5p_sandwich_mat_list[[3]])),
                              sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                              levels = c("SL", "2i", "mTORi"))
)

g1.1 <- ggplot(pol2s5p_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 21, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 221, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 1.5) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Pol II-S5p)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(21, 221), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "top")

g1.2 <- ggplot(ttseq_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 21, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 221, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 1.5) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Nascent RNA)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(21, 221), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "none")

g1.3 <- ggplot(TT_Pol2s5p_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 21, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 221, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 1.5) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Nascent RNA / Pol II-S5p)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(21, 221), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "none")

ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, ncol = 1, heights = c(4.7, 4, 4)),
       filename = "Fig3_TTseq_Pol2S5p_coverage.png", path = "figs",
       device = "png", width = 3.5, height = 10)

# PolII pausing index, start-seq TSS interval / (+500, +1500) gene body ----------------------
tss.gr <- importRanges("../data/tss.mm9.gff3")
tss.attr.gr <- importRanges("../data/tss.attributes.gff3")

library(EnsDb.Mmusculus.v79)
res <- biomaRt::select(EnsDb.Mmusculus.v79,
                       keys = as.character(tss.gr$gene_name),
                       keytype = "GENENAME",
                       columns = "GENEID")
tss.gr$gene_id <- res$GENEID[match(tss.gr$gene_name, res$GENENAME)]
tss.attr.gr$gene_id <- res$GENEID[match(tss.attr.gr$gene_name, res$GENENAME)]

Start_peaks.gr <- reduce(importRanges("../data/Start_seq_peaks.gtf"), min.gapwidth=50L)

gene_bodies.gr <- flank(promoters(tss.gr, upstream = 0, downstream = 1000),
                        width = 2000, start = F)

start_hits <- findOverlaps(Start_peaks.gr, tss.gr+100, ignore.strand = F) # adjust TSS to the nearest Start-seq peak
pause_sites.gr <- GRanges()
for(i in seq_along(tss.gr))
{
  tmp_hits <- queryHits(start_hits)[subjectHits(start_hits) == i]
  if (length(tmp_hits) == 1) {
    pause_sites.gr <- c(pause_sites.gr, Start_peaks.gr[tmp_hits])
  } else if (length(tmp_hits) > 1) {
    keep.idx <- nearest(tss.gr[i], Start_peaks.gr[tmp_hits], ignore.strand = F)
    pause_sites.gr <- c(pause_sites.gr, Start_peaks.gr[tmp_hits][keep.idx] )
  } else {
    pause_sites.gr <- c(pause_sites.gr, promoters(tss.gr[i], upstream = 0, downstream = 50))
  }
}
names(pause_sites.gr) <- names(tss.gr)
saveRDS(pause_sites.gr, "../fig4/data/tss_pause_sites.RData")

if (F) {
  # MINUTE ChIP data
  pausing_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1.*Pol5p.*ALL.*bam$", full.names = T),
                           intervals = pause_sites.gr, stranded = F, paired.end = "ignore")
  pausing_mat <- pausing_mat / width(pause_sites.gr) # convert to reads density

  gene_body_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1.*Pol5p.*ALL.*bam$", full.names = T),
                             intervals = gene_bodies.gr, stranded = F, paired.end = "ignore")
  gene_body_mat <- gene_body_mat / width(gene_bodies.gr) # convert to reads density
  
  pausing_index_mat <- pausing_mat / gene_body_mat
  
  pausing_index_mat <- pausing_index_mat[!is.na(tss.gr$gene_id), ]
  rownames(pausing_index_mat) <- tss.gr$gene_id[!is.na(tss.gr$gene_id)]
  
  pausing_index_mat <- pausing_index_mat[, c(3, 1, 4)]
  colnames(pausing_index_mat) <- c("SL", "2i", "mTORi")
  pausing_index_mat <- pausing_index_mat[!apply(pausing_index_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]
  pausing_index_mat <- log10(pausing_index_mat)
}

# est elongation speed heatmap  ---------------------------------------------------------
est_speed_SL_mat <- TT_Pol2s5p_sandwich_mat_list[[1]]

res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(rownames(est_speed_SL_mat)),
                       keytype = "ENTREZID",
                       columns = c("GENEID", "GENEBIOTYPE"))
rownames(est_speed_SL_mat) <- res$GENEID[match(rownames(est_speed_SL_mat), res$ENTREZID)]
est_speed_SL_mat <- est_speed_SL_mat[rownames(est_speed_SL_mat) %in% names(gene.gr), ]

est_speed_SL_mat_smooth <-
  apply(est_speed_SL_mat, 1, function(one)
    smooth.spline(x = 1:260, y = one, df = 15)$y) %>% t()

set.seed(1)
k <- kmeans(est_speed_SL_mat_smooth, 3)

speed_cls <- factor(k$cluster,
                    levels = order(aggregate(
                      rowMedians(est_speed_SL_mat_smooth, na.rm = T),
                      list(k$cluster),
                      "mean"
                    )[, 2])) %>% as.numeric() # Slow 4258 Medium 4106 Fast 2022 

# speed class ~ speed profile

speed_class_data <- data.frame(
  position = rep(seq_len(n_pos), 3),
  value = c(
    colMeans(est_speed_SL_mat[speed_cls == 1, ]),
    colMeans(est_speed_SL_mat[speed_cls == 2, ]),
    colMeans(est_speed_SL_mat[speed_cls == 3, ])
  ),
  sample = factor(c(
    rep("Slow", n_pos), rep("Medium", n_pos), rep("Fast", n_pos)
  ),
  levels = c("Slow", "Medium", "Fast"))
)

tss_speed_height <- speed_class_data[speed_class_data$position == 21, ]

ggplot(speed_class_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 21, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 220, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_n[c(1,10,6)]) +
  xlab("") +
  ylab("log Est. speed (a.u.)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(21, 221), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "none")

ggsave(filename = "Fig3_Est_speed_class_coverage.png", path = "figs",
       device = "png", width = 3.5, height = 3.5)


# speed class ~ pausing index
pausing_index_mat2 <- pausing_index_mat[rownames(est_speed_SL_mat), ]
data.frame(cls = as.factor(speed_cls),
           pi = log10(pausing_index_mat2[, 3])) %>%
  ggplot(aes(x = cls, y = pi, fill = as.factor(speed_cls), alpha = 1)) +
  geom_hline(yintercept = median(log10(pausing_index_mat2[speed_cls == 3, 3]), na.rm = T), 
             color = "grey50", lty = 2) +
  geom_violin() +
  geom_boxplot(outlier.color = NA, width = 0.3, size = 0.8) +
  ylim(c(-0.5, 3.3)) +
  xlab("") + ylab(expression('log'[10] * " Pausing index")) +
  scale_x_discrete(breaks = 1:3,
                   label = c("Slow", "Medium", "Fast")) +
  scale_fill_manual(values = colors_n[c(1,10,6)]) +
  theme_setting + 
  theme(legend.position = 'none')

ggsave(filename = "Fig3_Pausing_index_Est_speed_class.png", path = "figs",
       device = "png", width = 3.5, height = 3.5)

# heatmap
colorset = colorRampPalette(c("yellow", "white", "deepskyblue3", "black"))(100) #%>% rev()

plot_heatmap <- function(est_speed_SL_mat, cls = 1, .main = "") {
  mat <- est_speed_SL_mat[speed_cls == cls, ]
  image(t(mat[order(rowMedians(mat, na.rm = T)), ]),
        col = colorset,
        xlab='', xaxt='n',
        ylab='', yaxt='n',
        lwd=2, main = .main, cex = 5
  )
  box(col='black', lwd=2)
  axis(1, at=c(21, 221) / 260, labels=c("TSS","TTS"), cex = 6 )
  axis(2, at=c(nrow(mat), 1) / nrow(mat), labels=c(1, nrow(mat)), cex=1 )
}

if (T) {
  png(filename = "../figS3/figs/FigS3_Est_speed_class_heatmap.png", width = 800, height = 400)
  
  layout(mat = matrix(1:4, 1),
         widths = c(3,3,3,0.5),
         heights = c(2,2,2,2)+0.25, TRUE)
  layout.show(n = 4)
  
  par(mar=c(1,2,1,1))
  
  plot_heatmap(est_speed_SL_mat, 1, "Slow")
  plot_heatmap(est_speed_SL_mat, 2, "Medium")
  plot_heatmap(est_speed_SL_mat, 3, "Fast")
  
  # scales
  par(mar=c(1,0,1,2))
  image(x = 1, y = seq(min(est_speed_SL_mat), max(est_speed_SL_mat), length.out=100), 
        matrix(seq(min(est_speed_SL_mat), max(est_speed_SL_mat), length.out=100), 1, 100),
        col = colorset,
        xlab='', ylab='',
        main='', xaxt='n', yaxt='n',
        lwd=1, axes=TRUE)
  axis(4)
  box(col='black', lwd=1.5)
  
  dev.off()
}

# plot Pol2s5p density tss ~ gene body --------------------------------------------------
if (F) {
  dat <- data.frame(TSS = log10(pausing_mat[,3]),
                    Gene_body = log10(gene_body_mat[,3])) %>%
    dplyr::filter(is.finite(TSS) & is.finite(Gene_body))
  g0.1 <- dat %>%
    ggplot(aes(x = TSS, y = Gene_body)) +
    geom_hex(binwidth=c(0.05, 0.08)) +
    geom_text(x = 0.6, y = 0, label = paste("r =", round(cor(dat$TSS, dat$Gene_body), 3) )) +
    geom_text(x = 0.6, y = -0.2, label = paste("n =", nrow(dat) )) +
    theme_setting +
    scale_fill_viridis_c(option = "D", direction = -1, end = 0.93) +
    scale_x_continuous(name="\nTSS density", breaks=c((-2):1), labels=10^((-2):1), limits = c(-2, 1)) +
    scale_y_continuous(name="Gene body density", breaks=c((-3):0), labels=10^((-3):0), limits = c(-3, 0)) +
    labs(fill = "Genes", title = "Pol II-S5p SL")+
    theme(legend.position = "none")
  
  dat <- data.frame(TSS = log10(pausing_mat[,1]),
                    Gene_body = log10(gene_body_mat[,1])) %>%
    dplyr::filter(is.finite(TSS) & is.finite(Gene_body))
  g0.2 <- dat %>%
    ggplot(aes(x = TSS, y = Gene_body)) +
    geom_hex(binwidth=c(0.05, 0.08)) +
    geom_text(x = 0.6, y = 0, label = paste("r =", round(cor(dat$TSS, dat$Gene_body), 3) )) +
    geom_text(x = 0.6, y = -0.2, label = paste("n =", nrow(dat) )) +
    theme_setting +
    scale_fill_viridis_c(option = "D", direction = -1, end = 0.93) +
    scale_x_continuous(name="\nTSS density", breaks=c((-2):1), labels=10^((-2):1), limits = c(-2, 1)) +
    scale_y_continuous(name="Gene body density", breaks=c((-3):0), labels=10^((-3):0), limits = c(-3, 0)) +
    labs(fill = "Genes", title = "Pol II-S5p 2i")+
    theme(legend.position = "none")
  
  dat <- data.frame(TSS = log10(pausing_mat[,4]),
                    Gene_body = log10(gene_body_mat[,4])) %>%
    dplyr::filter(is.finite(TSS) & is.finite(Gene_body))
  g0.3 <- dat %>%
    ggplot(aes(x = TSS, y = Gene_body)) +
    geom_hex(binwidth=c(0.05, 0.08)) +
    geom_text(x = 0.6, y = 0, label = paste("r =", round(cor(dat$TSS, dat$Gene_body), 3) )) +
    geom_text(x = 0.6, y = -0.2, label = paste("n =", nrow(dat) )) +
    theme_setting +
    scale_fill_viridis_c(option = "D", direction = -1, end = 0.93) +
    scale_x_continuous(name="\nTSS density", breaks=c((-2):1), labels=10^((-2):1), limits = c(-2, 1)) +
    scale_y_continuous(name="Gene body density", breaks=c((-3):0), labels=10^((-3):0), limits = c(-3, 0)) +
    labs(fill = "Genes", title = "Pol II-S5p mTORi")+
    theme(legend.position = "none")
  
  ggsave(plot = grid.arrange(g0.1, g0.2, g0.3, nrow = 1),
         filename = "FigS3_Pol2S5p_TSS_gene_body.png", path = "../figS3/figs",
         device = "png", width = 10, height = 3.5)
}

# plot pause index ~ RNA synthesis --------------------------------------------------
pausing_index_mat <- readRDS("../fig3/data/pausing_index_SL_2i.RData")
pausing_index_mat <- pausing_index_mat[order(pausing_index_mat[,1], decreasing = T), c(3,1,4)]

# add Tx measurements
sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$gene_id %in% rownames(pausing_index_mat), ]
sample_ord <- match(rownames(pausing_index_mat), unique(sample_Tx_counts$gene_id))

pausing_production_mat <- pausing_index_mat
for ( i in c("SL", "2i_2d", "mTORi_2d")) {
  pausing_production_mat <- 
    cbind(pausing_production_mat, # append nascent RNA production speed = copy * labeled rate
          rowSums(log10(sample_Tx_counts[sample_Tx_counts$Sample == i, 
                                         c("Copy", "Labeled_rate")][sample_ord, ])) - log10(5))
}
colnames(pausing_production_mat) <- c("SL", "2i_2d", "mTORi_2d",
                                      paste0("mu_", c("SL", "2i_2d", "mTORi_2d")))
pausing_production_mat <- 
  pausing_production_mat[!apply(pausing_production_mat, 1,
                                function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# negative correlation between PI and gene expression:
# https://doi.org/10.1093/nar/gkx1225
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0984-2

g1 <- data.frame(pause_index = pausing_production_mat[, 1], # n = 4492
                 RNA_production = pausing_production_mat[, 4]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(bins= 50) +
  geom_hline(yintercept = median(pausing_production_mat[, 1]),
             lty = 2,
             alpha = 0.5) +
  geom_text(x = -0.25, 
            y = 3, 
            label = paste("r =", round(cor(pausing_production_mat[, c(1,4)])[1,2], 2) )) +
  geom_text(x = -0.2,
            y = median(pausing_production_mat[, 1]) + 0.1, 
            label = "median") +
  geom_text(x = -0.2,
            y = median(pausing_production_mat[, 1]) - 0.1,
            label = round(10^median(pausing_production_mat[, 1]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (copy/min)",
                     breaks=c(-3, -2, -1, 0), 
                     labels=c(0.001, 0.01, 0.1, 1), 
                     limits = c(-3, 0.2)) +
  scale_y_continuous(name="Pausing Index", 
                     breaks=0:3, 
                     labels=10^(0:3),
                     limits = c(-0.5, 3)) +
  labs(fill = "Genes", title = "SL")+
  theme(legend.position = "none")

g2 <- data.frame(pause_index = pausing_production_mat[, 2],
           RNA_production = pausing_production_mat[, 5]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(binwidth=c(0.05, 0.08)) +
  geom_hline(yintercept = median(pausing_production_mat[, 2]), lty = 2, alpha = 0.5) +
  geom_text(x = -0.25, 
            y = 3, 
            label = paste("r =", round(cor(pausing_production_mat[, c(2,5)])[1,2], 2) )) +
  geom_text(x = -0.2, 
            y = median(pausing_production_mat[, 2]) + 0.1, 
            label = "median") +
  geom_text(x = -0.2, 
            y = median(pausing_production_mat[, 2]) - 0.1, 
            label = round(10^median(pausing_production_mat[, 2]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (copy/min)",
                     breaks=c(-3, -2, -1, 0), 
                     labels=c(0.001, 0.01, 0.1, 1), 
                     limits = c(-3, 0.2)) +
  scale_y_continuous(name="Pausing Index", 
                     breaks=0:3, labels=10^(0:3),
                     limits = c(-0.5, 3)) +
  labs(fill = "Genes", title = "2i")+
  theme(legend.position = "none")

g3 <- data.frame(pause_index = pausing_production_mat[, 3],
           RNA_production = pausing_production_mat[, 6]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(binwidth=c(0.05, 0.08)) +
  geom_hline(yintercept = median(pausing_production_mat[, 3]),
             lty = 2, 
             alpha = 0.5) +
  geom_text(x = -0.25,
            y = 3, 
            label = paste("r =", round(cor(pausing_production_mat[, c(3,6)])[1,2], 2) )) +
  geom_text(x = -0.2, 
            y = median(pausing_production_mat[, 3]) + 0.1,
            label = "median") +
  geom_text(x = -0.2, 
            y = median(pausing_production_mat[, 3]) - 0.1, 
            label = round(10^median(pausing_production_mat[, 3]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (copy/min)",
                     breaks=c(-3, -2, -1, 0),
                     labels=c(0.001, 0.01, 0.1, 1), 
                     limits = c(-3, 0.2)) +
  scale_y_continuous(name="Pausing Index",
                     breaks=0:3,
                     labels=10^(0:3), 
                     limits = c(-0.5, 3)) +
  labs(title = "mTORi") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "FigS3_Pause_index_tx_correlation.png",
       path = "../figS3/figs",
       device = "png", width = 10, height = 3.5)

# plot gene body Pol II ~ RNA synthesis --------------------------------------------------

g4 <- ggplot(gene_body_production_mat, aes(x = Tx_RPK_SL, y = Pol2_SL)) + # n = 5890
  # geom_hex(bins= 50) +
  geom_point(aes(color = with(gene_body_production_mat, get_dens(Tx_RPK_SL, Pol2_SL)))) + 
  geom_hline(yintercept = median(gene_body_production_mat[, 1]), lty = 2, alpha = 0.5) +
  annotate(geom = "text", x = 2.5, y = 2.8, hjust = "left", vjust = "top",
            label = paste0(paste("r = ", round(cor(gene_body_production_mat[, c(1,4)])[1,2], 2)),
                           "\nn = ", nrow(gene_body_production_mat))) +
  geom_text(x = -0.5, 
            y = median(gene_body_production_mat[, 1]) + 0.1, 
            label = "median") +
  geom_text(x = -0.5,
            y = median(gene_body_production_mat[, 1]) - 0.1,
            label = round(10^median(gene_body_production_mat[, 1]), 1)) +
  theme_setting +
  scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (RPK)", 
                     breaks=(0:3),
                     labels=10^(0:3), 
                     limits = c(-1, 3.8)) +
  scale_y_continuous(name="Pol II S5p gene body (RPK)",
                     breaks=0:2, labels=10^(0:2),
                     limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "SL")+
  theme(legend.position = "none")

g5 <- ggplot(gene_body_production_mat, aes(x = Tx_RPK_2i, y = Pol2_2i)) + # n = 5890
  geom_point(aes(color = with(gene_body_production_mat, get_dens(Tx_RPK_2i, Pol2_2i)))) + 
  geom_hline(yintercept = median(gene_body_production_mat[, 2]), lty = 2, alpha = 0.5) +
  annotate(geom = "text", x = 2.5, y = 2.8, hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(cor(gene_body_production_mat[, c(2,5)])[1,2], 2)),
                          "\nn = ", nrow(gene_body_production_mat))) +
  geom_text(x = -0.5,
            y = median(gene_body_production_mat[, 2]) + 0.1,
            label = "median") +
  geom_text(x = -0.5, 
            y = median(gene_body_production_mat[, 2]) - 0.1, 
            label = round(10^median(gene_body_production_mat[, 2]), 1)) +
  theme_setting +
  scale_color_viridis_c(option = "A", 
                       direction = -1,
                       begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (RPK)", 
                     breaks=(0:3),
                     labels=10^(0:3), 
                     limits = c(-1, 3.8)) +
  scale_y_continuous(name="Pol II S5p gene body (RPK)", 
                     breaks=0:2, 
                     labels=10^(0:2), 
                     limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "2i")+
  theme(legend.position = "none")

g6 <- ggplot(gene_body_production_mat, 
             aes(x = Tx_RPK_mTORi, y = Pol2_mTORi)) + # n = 5890
  geom_point(aes(color = with(gene_body_production_mat, get_dens(Tx_RPK_mTORi, Pol2_mTORi)))) + 
  geom_hline(yintercept = median(gene_body_production_mat[, 3]), lty = 2, alpha = 0.5) +
  annotate(geom = "text", x = 2.5, y = 2.8, hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(cor(gene_body_production_mat[, c(3,6)])[1,2], 2)),
                          "\nn = ", nrow(gene_body_production_mat))) +
  geom_text(x = -0.5, 
            y = median(gene_body_production_mat[, 3]) + 0.1, 
            label = "median") +
  geom_text(x = -0.5,
            y = median(gene_body_production_mat[, 3]) - 0.1, 
            label = round(10^median(gene_body_production_mat[, 3]), 1)) +
  theme_setting +
  scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nTranscription (RPK)", 
                     breaks=(0:3),
                     labels=10^(0:3), 
                     limits = c(-1, 3.8)) +
  scale_y_continuous(name="Pol II S5p gene body (RPK)", 
                     breaks=0:2, labels=10^(0:2), 
                     limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "mTORi")+
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g4, g5, g6, ncol = 1),
       filename = "Fig3_Pol2S5p_genebody_tx.png", path = "figs",
       device = "png", width = 3.5, height = 10)

# plot RNA synthesis explained ~ features --------------------------------------------------
if (F) {
  genes.intercect <- intersect.Vector(rownames(readRDS("../fig5/data/elongation_speed_table.RData")), 
                                      intersect.Vector(rownames(pausing_index_mat), rownames(gene_body_production_mat)))
  
  dat <- cbind(pausing_production_mat[genes.intercect, 1:3],
               gene_body_production_mat[genes.intercect, ],
               DHS = log(as.numeric(tss.attr.gr$DHS_density[match(genes.intercect, tss.attr.gr$gene_id)]))
  )
  
  dat$CpG <- interval_pattern_hits(intervals = promoters(gene.gr[match(genes.intercect, gene.gr$gene_id)]),
                                   pattern = "CG", 
                                   which_genome = "mm9") %>% log1p
  dat$TATA <- interval_pattern_hits(intervals = promoters(gene.gr[match(genes.intercect, gene.gr$gene_id)]),
                                    pattern = "TATA", 
                                    which_genome = "mm9") %>% log1p
  
  R2_data <- data.frame(R2 = c(multi_variance_explained(trim_quantile(dat[, c("SL", "Pol2_SL", "CpG", "TATA", "DHS")]), trim_quantile(dat$tx_SL)),
                               multi_variance_explained(trim_quantile(dat[, c("2i", "Pol2_2i", "CpG", "TATA", "DHS")]), trim_quantile(dat$tx_2i)),
                               multi_variance_explained(trim_quantile(dat[, c("mTORi", "Pol2_mTORi", "CpG", "TATA", "DHS")]), trim_quantile(dat$tx_mTORi))),
                        Sample = rep(c("SL", "2i", "mTORi"), each = 5),
                        Feature = rep(c("Pausing index", "Pol2S5p GB", "CpG", "TATA", "DHS"), 3))
  R2_data$Sample <- factor(R2_data$Sample, levels = c("SL", "2i", "mTORi"))
  R2_data$Feature <- factor(R2_data$Feature, c("CpG", "TATA", "Pausing index", "DHS", "Pol2S5p GB"))
  
  ggplot(R2_data, aes(x = Sample, y = R2, fill = Feature)) +
    geom_bar(position="dodge", stat = "identity", width = 0.9) +
    xlab("R-squared composition") +
    ylab("RNA synthesis explained") +
    scale_fill_manual(values = colors_9[c(2,8,1,5,4)]) +
    theme_setting +
    theme(legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  ggsave(filename = "Fig3_RNA_synthesis_explained.png", path = "figs",
         device = "png", width = 4, height = 3)
}

# plot estimated speed ~ measured speed --------------------------------------------------
elongation_speed_table <- read.table("../fig3/data/elongation_speed_table.txt")

genes.intercect <- intersect.Vector(rownames(elongation_speed_table), 
                                    rownames(gene_body_production_mat))
dat_speed <- with(gene_body_production_mat[genes.intercect, ],
                  data.frame(# estimated speed
                             v_hat_SL = (Tx_RPK_SL - Pol2_SL),
                             v_hat_2i = (Tx_RPK_2i - Pol2_2i),
                             v_hat_mTORi = (Tx_RPK_mTORi - Pol2_mTORi), 
                             v_hat_all = (LRNA_RPK - Pol2_mTORi),
                             # est. synthesis, mu, and synthesis mu == LRNA_RPK
                             mu_hat = log10(elongation_speed_table[genes.intercect, "speed"]) + Pol2_SL,
                             mu = pausing_production_mat[genes.intercect, "mu_SL"],
                             # measured speed
                             speed = log10(elongation_speed_table[genes.intercect, "speed"]),
                             Pol2_den = (Pol2_SL)
                             )
                  )
rownames(dat_speed) <- genes.intercect
dat_speed <- dat_speed %>% dplyr::filter(speed > -1)


if (T) {
  # validate measured speed from Gro-seq flv inhibition with provided table
  genes.intercect2 <- match(paper_est_speed$gene_id, rownames(elongation_speed_table))
  
  dat_speed_cmp2 <- data.frame(speed_paper = paper_est_speed[, "speed"],
                               speed_measured = elongation_speed_table[genes.intercect2, "speed"],
                               Estimation_at = paste(paper_est_speed[, "Time"], "min"),
                               gene_id = paper_est_speed$gene_id) 
  dat_speed_cmp2 <- dat_speed_cmp2[!is.na(dat_speed_cmp2$speed_measured), ]
  dat_speed_cmp2$Estimation_at <- factor(dat_speed_cmp2$Estimation_at)
  
  ggplot(dat_speed_cmp2, aes(x = (speed_paper), y = (speed_measured))) +
    geom_abline(slope = 1, intercept = 0, color = "grey50") +
    geom_point(aes(color = Estimation_at), size = 1) + 
    scale_x_continuous(name="Speed Jonker et al. (kb/min)",
                       breaks= c(0, 1, 2, 3),
                       labels= c(0, 1, 2, 3), limits = c(0, 3.5)) +
    scale_y_continuous(name="Measured speed (kb/min)",
                       breaks= c(0, 1, 2, 3),
                       labels= c(0, 1, 2, 3)) +
    theme_setting 
  ggsave(filename = paste0("FigS4_measured_speed_vs_paper_speed.png"), 
         path = "../figS4/figs/",
         device = "png", width = 6, height = 4.5)
}

# plot speed, TX and Pol II interplays
g7 <- ggplot(dat_speed, aes(x = v_hat_SL, y = speed)) +
  geom_point(aes(color = get_dens(v_hat_SL, speed)), 
             size = 0.8) +
  scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
  annotate(geom = "text", x = 1.3, y = 0.85, hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(cor(dat_speed$v_hat_SL, dat_speed$speed), 3)),
                          "\nn = ", nrow(dat_speed))) +
  scale_x_continuous(name="\nEstimated velocity (a.u.)\n",
                     limits = c(-0.2, 2)) +
  scale_y_continuous(name="\nMeasured velocity (kb/min)",
                     breaks=c(-0.7, 0, 0.7),
                     labels=round(10^c(-0.7, 0, 0.7), 1),
                     limits = c(-1, 0.9)) +
  theme_setting +
  theme(legend.position = "none") 

# plot estimated speed ~ samples ------------------------------------------------------------
dat_speed_cmp = data.frame(Sample = rep(c("SL", "2i", "mTORi"), 
                                        each = nrow(dat_speed)),
                           Speed_hat = c(as.matrix(dat_speed[, 1:3])))
dat_speed_cmp$Sample <- factor(dat_speed_cmp$Sample, levels = rev(c("SL", "2i", "mTORi")))

g8 <- ggplot(dat_speed_cmp,
             aes(x = Speed_hat, y = Sample, fill = Sample)) + 
  ggridges::geom_density_ridges(rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2, lty = 2) +
  theme_setting +
  xlab("\nEstimated velocity (a.u.)\n") +
  ylab("") +
  scale_fill_manual(values = colors_20[rev(c(13, 2, 7))]) +
  theme(legend.position = "none")
  
# plot RNA synthesis ~ estimated synthesis --------------------------------------------------
if (F) { # previous version
  g10 <- ggplot(dat_speed[complete.cases(dat_speed[, "mu"]), ], 
                aes( x = mu_hat, y = mu)) +
    geom_point(aes(color = get_dens(mu, mu_hat)), size = 0.8) + 
    scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
    annotate(geom = "text", x = 1.4, y = 0.1, hjust = "left", vjust = "top",
             label = paste0(paste("r = ", 
                                  round(single_variance_explained(dat_speed$mu_hat, dat_speed$mu, T), 3),
                                  "\nn = ", nrow(dat_speed))) ) +
    xlab("\nEst. RNA synthesis (a.u.)\n") +
    xlim(c(0, 2)) +
    scale_y_continuous(name="\nRNA synthesis (copy/min)",
                       breaks=(-2:0), labels=10^(-2:0), limits = c(-2, 0.1)) +
    theme_setting +
    theme(legend.position = "none")
}

if (T) {  # mu ~ mu_hat
  # consider all measured velocity by any number of time points (Flv GRO-seq)
  # the estimated RNA synthesis still remains a weak correlation with actual RNA synthesis rates 
  # but RNA synthesis only marginally correlates with the measured velocity by more than 2 time points, which are the genes with longer lengths and better GRO-seq coverages
  elongation_speed_table <- read.table("../fig3/data/elongation_speed_table.txt")
  
  sample_Tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined.RData")
  sample_Tx_counts_Rates <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$gene_id %in% rownames(gene_body_production_mat) &
                                                     sample_Tx_counts_Rates$Sample == "SL", ]
  
  dat_mu_hat <- data.frame(mu_hat = gene_body_production_mat$Pol2_SL + 
                             log10(elongation_speed_table[match(rownames(gene_body_production_mat), rownames(elongation_speed_table)), 1]),
                           mu = sample_Tx_counts_Rates[match(rownames(gene_body_production_mat), sample_Tx_counts_Rates$gene_id), c("Copy", "Labeled_rate")] %>% 
                             log10() %>% rowSums() %>% "-"(log10(5))
  ) %>% dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>% 
    dplyr::filter(mu_hat > 0.2 & mu_hat < 2 & mu < 2 )
  
  g10 <- ggplot(dat_mu_hat, aes(x = mu_hat, y = mu, color = get_dens(mu_hat, mu))) +
    geom_point() +
    annotate(geom = "text", x = Inf, y = Inf, hjust = 1, vjust = 1.2,
             label = paste0(paste("r = ", 
                                  with(dat_mu_hat, 
                                       single_variance_explained(mu, mu_hat, T)) %>% round(3),
                                  "\nn = ", nrow(dat_mu_hat))) ) +
    scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
    scale_y_continuous(name = "\nRNA synthesis (copy/min)",
                       breaks = c(-2:1), labels = 10^(-2:1), 
                       limits = c(-2, 0.5)) +
    scale_x_continuous(name = "\nEst. RNA synthesis (a.u.)\n", 
                       breaks = c(0:2), labels = 10^(0:2),
                       limits = c(0, 2.2)) +
    theme_setting +
    theme(legend.position = "none")
}

# plot Pol II ~ speed -------------------------------------------------------
g11 <- ggplot(dat_speed, aes(x = Pol2_den, y = speed)) +
  geom_point(aes(color = get_dens(Pol2_den, speed)), size = 0.8) + 
  scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
  annotate(geom = "text", x = 1.4, y = 0.8, hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(cor(dat_speed$speed, dat_speed$Pol2_den), 3)),
                          "\nn = ", nrow(dat_speed))) +
  scale_x_continuous(name="\nPol II S5p density (RPK)\n",
                     breaks=c(0, 0.5, 1, 1.5, 2),
                     labels=round(10^c(0, 0.5, 1, 1.5, 2), 0),
                     limits = c(0, 2)) +
  scale_y_continuous(name="\nMeasured velocity (kb/min)",
                     breaks=c(-0.7, 0, 0.7),
                     labels=round(10^c(-0.7, 0, 0.7), 1),
                     limits = c(-1, 0.8)) +
  theme_setting +
  theme(legend.position = "none") 

# plot RNA synthesis ~ measured speed --------------------------------------------------
if (T) { # previous version
  g12 <- ggplot(dat_speed[complete.cases(dat_speed[, "mu"]), ], 
                aes(x = mu, y = speed)) +
    geom_point(aes(color = get_dens(mu, speed)), size = 0.8) + 
    scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
    annotate(geom = "text", x = -0.5, y = 0.8, hjust = "left", vjust = "top",
             label = paste0(paste("r = ", 
                                  round(single_variance_explained(dat_speed$speed, dat_speed$mu, T), 3)),
                            "\nn = ", nrow(dat_speed))) +
    scale_x_continuous(name="\nRNA synthesis (copy/min)\n",
                       breaks=(-2:0), labels=10^(-2:0), limits = c(-2, 0.1)) +
    scale_y_continuous(name="\nMeasured velocity (kb/min)",
                       breaks=c(-0.7, 0, 0.7),
                       labels=round(10^c(-0.7, 0, 0.7), 1),
                       limits = c(-0.8, 0.8)) +
    theme_setting +
    theme(legend.position = "none") 
}

if (F) { # revision, 2021 Jul 26
  # velocity ~ RNA copy
  elongation_speed_table <- read.table("../fig3/data/elongation_speed_table.txt")
  # elongation_speed_table <- elongation_speed_table[elongation_speed_table$time_points > 1, ]
  
  sample_Tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined.RData")
  sample_Tx_counts_Rates <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$gene_id %in% rownames(elongation_speed_table) &
                                                     sample_Tx_counts_Rates$Sample == "SL", ]
  sample_Tx_counts_Rates$Speed <- elongation_speed_table[match(sample_Tx_counts_Rates$gene_id, rownames(elongation_speed_table)), 1]
  sample_Tx_counts_Rates <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$Copy > 0, ]
  sample_Tx_counts_Rates <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$Labeled_rate > 0, ]
  sample_Tx_counts_Rates$Copy <- log10(sample_Tx_counts_Rates$Copy)
  sample_Tx_counts_Rates$Speed <- log10(sample_Tx_counts_Rates$Speed)
  sample_Tx_counts_Rates$mu <- sample_Tx_counts_Rates$Copy + log10(sample_Tx_counts_Rates$Labeled_rate) - log10(5)
  sample_Tx_counts_Rates <- sample_Tx_counts_Rates %>% 
    dplyr::filter(complete.cases(.)) %>% dplyr::filter(Copy > -3)
  
  ggplot(sample_Tx_counts_Rates, aes(x = Speed, y = Copy, color = get_dens(Speed, Copy))) +
    geom_point() +
    annotate(geom = "text", x = Inf, y = -Inf, hjust = 1.5, vjust = -1,
             label = paste0(paste("r = ", 
                                  with(sample_Tx_counts_Rates, 
                                       single_variance_explained(Copy, Speed, T)) %>% round(3),
                                  "\nn = ", nrow(sample_Tx_counts_Rates))) ) +
    scale_color_viridis(option = "C", direction = -1, alpha = 0.8) +
    scale_x_continuous(name = "\nMeasured velocity (kb/min)\n",
                       breaks = c(-0.7, 0, 0.7), labels = c(0.2, 1, 5),
                       limits = c(-1, 0.7)) +
    scale_y_continuous(name = "\nTotal mRNA Copy", breaks = c(-1:2), labels = 10^(-1:2)) +
    theme_setting +
    theme(legend.position = "none")
  
  ggsave(filename = "FigS3_scatter_copy_measured_speed.png",
         path = "../figS3/figs",
         device = "png", width = 4, height = 4)
}

ggsave(plot = cowplot::plot_grid(g11, g7, g12, g10, g8, ncol = 2, align = "v"),
       filename = "Fig3_TTseq_Pol2S5p_speed_interplays2.png", path = "../fig3/figs",
       device = "png", width = 7, height = 10)


# plot speed changes in 2i / mTORi ----------------------------------------------------
dat_speed_change <- with(gene_body_production_mat,
                  data.frame(# estimated speed
                    v_hat = ((Tx_RPK_SL - Pol2_SL) + (Tx_RPK_2i - Pol2_2i) + (Tx_RPK_mTORi - Pol2_mTORi)) / 3,
                    log10FC_2i = (Tx_RPK_2i - Pol2_2i) - (Tx_RPK_SL - Pol2_SL),
                    log10FC_mTORi = (Tx_RPK_mTORi - Pol2_mTORi) - (Tx_RPK_SL - Pol2_SL),
                    log10FC_2i_mTORi = (Tx_RPK_2i - Pol2_2i) - (Tx_RPK_mTORi - Pol2_mTORi) 
                  )
)
rownames(dat_speed_change) <- rownames(gene_body_production_mat)
saveRDS(dat_speed_change, "../fig3/data/dat_speed_change.RData")

dat_speed_change <- dat_speed_change %>% 
  dplyr::filter(complete.cases(.) & v_hat > (-1) & !is.infinite(v_hat) ) %>%
  dplyr::mutate("density_2i" = get_dens(v_hat, log10FC_2i, n.grid = 200)) %>%
  dplyr::mutate("density_mTORi" = get_dens(v_hat, log10FC_mTORi, n.grid = 200)) 

g13 <- ggplot(dat_speed_change, aes(x = v_hat, y = log10FC_2i, color = density_2i)) +
  geom_hline(yintercept = 0,
             lwd = 1, 
             color = add.alpha("blue", 0.5)) +
  geom_point(size = 1) + 
  geom_smooth(method = "gam",
              lwd = 0.5, 
              color = "red2",
              se = F) +
  scale_color_gradient(low = add.alpha("grey50", 0.5),
                       high = "orange") +
  ylim(c(-3, 1)) +
  xlab("log10 Est. speed (a.u.)") +
  ylab("log10FC") +
  ggtitle("2i speed change") +
  theme_setting +
  theme(legend.position = "none")

g14 <- ggplot(dat_speed_change, aes(x = v_hat, y = log10FC_mTORi, color = density_mTORi)) +
  geom_hline(yintercept = 0,
             lwd = 1, 
             color = add.alpha("blue", 0.5)) +
  geom_point(size = 1) + 
  geom_smooth(method = "gam",
              lwd = 0.5, 
              color = "red2",
              se = F) +
  scale_color_gradient(low = add.alpha("grey50", 0.5),
                       high = "green4") +
  ylim(c(-3, 1)) +
  xlab("log10 Est. speed (a.u.)") +
  ylab("log10FC") +
  ggtitle("mTORi speed change") +
  theme_setting +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g13, g14, ncol = 2),
       filename = "FigS4_MAplot_speed_change.png", 
       path = "../figS4/figs", device = "png", width = 12, height = 6)

# -------------------------------------------------------------------------------------- #
# SL2i velocity coverage
n_pos <- ncol(TT_Pol2s5p_sl2i_sandwich_mat[[1]])
TT_Pol2s5p_sl2i_data <- data.frame(position = rep(seq_len(n_pos), 2),
                              value = c(colMeans(TT_Pol2s5p_sl2i_sandwich_mat[[1]]),
                                        colMeans(TT_Pol2s5p_sl2i_sandwich_mat[[2]])),
                              sample = factor(c(rep("SL", n_pos), rep("SL2i", n_pos)),
                                              levels = c("SL", "SL2i"))
)

ggplot(TT_Pol2s5p_sl2i_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 21, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 221, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 1.5) +
  scale_color_manual(values = colors_20[c(13, 10)]) +
  xlab("") +
  ylab("log(Nascent RNA / Pol II-S5p)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(21, 221), labels = c("TSS", "TES")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "top")

ggsave(filename = "FigS3_SL_SL2i_TTseq_Pol2S5p_coverage.png", 
       path = "../figS3/figs", device = "png", width = 3.5, height = 3.5)

# -------------------------------------------------------------------------------------- #





