#-----------------------------------------------------------------------------------------------------
# This part estimate transcription elongation dynamics with TTseq and Pol II S5p coverage
# evaluate RNA Pol II pausing in different mouse ES pluripotent states
#
# Rui Shao, June 2020
#-----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

#-----------------------------------------------------------------------------------------------------
# estimated speed from TT-seq and Pol II S5p
# Define of parameters:
#   v_bar: average elongation speed
#   l: gene length
#   P_0: number of Pol II starting elongation per unit of time
#   P_1: total number of Pol II on gene body engaging on elongation
#   P_bar: number of Pol II per Kb
#   S_bar: RNA synthesis rate (~ TT-seq LRNA RPK or Copy * labeled rate)
#   v_hat: estimated mean speed 
#   m: minutes of elongation

# Equations 
# 1) P_0 ~ S_bar (~ links with a constant between ChIP and TT-seq)
# 2) P_1 = P_0 * l / v_bar
# 3) local Pol II speed v_i = P0 / d_i, where d_i is local Pol II RPK
# 4) v_bar = \sum_i^m v_i / m
# 5) P_bar = P_1 / l = P_0 / v_bar

#-----------------------------------------------------------------------------------------------------
# use txdb gene references
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(gene.gr$gene_id),
                       keytype = "ENTREZID",
                       columns = "GENEID")
gene.gr$gene_id <- res$GENEID[match(gene.gr$gene_id, res$ENTREZID)]
# keep only active genes, longer than 2kb
gene.gr <- gene.gr[width(gene.gr) > 2000 & width(gene.gr) < 1000000]
gene.gr <- gene.gr[seqnames(gene.gr) %in% paste0("chr", c(1:19, "X", "Y"))]
txRPK <- readRDS("../data/txRPK_SL_2i.RData")
gene.gr <- sort(gene.gr[na.omit(match(rownames(txRPK)[rowMeans(txRPK) > 0.2], gene.gr$gene_id))])
saveRDS(gene.gr, "../fig4/data/gene.gr.RData")

# TU.gr <- importRanges("../data/TU_anno/TU_filter+LRNA_SL_mm9.gtf")
# TU.gr$gene_id <- gsub("\\..", "\\1", TU.gr$gene_id)

tss.gr <- importRanges("../data/tss.mm9.gff3")

library(EnsDb.Mmusculus.v79)
res <- biomaRt::select(EnsDb.Mmusculus.v79,
                       keys = as.character(tss.gr$gene_name),
                       keytype = "GENENAME",
                       columns = "GENEID")
tss.gr$gene_id <- res$GENEID[match(tss.gr$gene_name, res$GENENAME)]

Start_peaks.gr <- reduce(importRanges("../data/Start_seq_peaks.gtf"), min.gapwidth=50L)

######################################################################################################
# PolII pausing index, start-seq TSS interval / (+500,+1500) gene body
######################################################################################################
gene_bodies.gr <- flank(promoters(tss.gr, upstream = 0, downstream = 1000), width = 2000, start = F)

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
names(pause_sites.gr) <- NULL

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
  
  # saveRDS(pausing_index_mat, "data/pausing_index_mattx.tss_SL_2i.RData")
}

# plot Pol2s5p density tss ~ gene body --------------------------------------------------
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
       filename = "Fig4_Pol2S5p_TSS_gene_body.png", path = "figs",
       device = "png", width = 10, height = 3.5)

# plot pause index ~ RNA synthesis --------------------------------------------------
pausing_index_mat <- pausing_index_mat[, c(3, 1, 4)]
colnames(pausing_index_mat) <- c("SL", "2i", "mTORi")
pausing_index_mat <- pausing_index_mat[!apply(pausing_index_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]
pausing_index_mat <- log10(pausing_index_mat)

pausing_index_mat <- pausing_index_mat[order(pausing_index_mat[,1], decreasing = T), ]

# add Tx measurements
sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$gene_id %in% rownames(pausing_index_mat), ]

sample_ord <- match(rownames(pausing_index_mat), unique(sample_Tx_counts$gene_id))

pausing_production_mat <- pausing_index_mat
for ( i in c("SL", "2i_2d", "mTORi_2d")) {
  pausing_production_mat <- cbind(pausing_production_mat, # append nascent RNA production speed = copy * labeled rate
                              rowSums(log10(sample_Tx_counts[sample_Tx_counts$Sample == i, c("Copy", "Labeled_rate")][sample_ord, ])) - log10(5))
}
pausing_production_mat <- pausing_production_mat[!apply(pausing_production_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# negative correlation between PI and gene expression:
# https://doi.org/10.1093/nar/gkx1225
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0984-2


g1 <- data.frame(pause_index = pausing_production_mat[, 1], # n = 4492
                 RNA_production = pausing_production_mat[, 4]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(bins= 50) +
  geom_hline(yintercept = median(pausing_production_mat[, 1]), lty = 2, alpha = 0.5) +
  geom_text(x = -0.25, y = 3, label = paste("r =", round(cor(pausing_production_mat[, c(1,4)])[1,2], 2) )) +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 1]) + 0.1, label = "median") +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 1]) - 0.1, label = round(10^median(pausing_production_mat[, 1]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3, -2, -1, 0), labels=c(0.001, 0.01, 0.1, 1), limits = c(-3, 0)) +
  scale_y_continuous(name="Pausing Index", breaks=0:3, labels=10^(0:3), limits = c(-0.5, 3)) +
  labs(fill = "Genes", title = "SL")+
  theme(legend.position = "none")

g2 <- data.frame(pause_index = pausing_production_mat[, 2],
           RNA_production = pausing_production_mat[, 5]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(binwidth=c(0.05, 0.08)) +
  geom_hline(yintercept = median(pausing_production_mat[, 2]), lty = 2, alpha = 0.5) +
  geom_text(x = -0.25, y = 3, label = paste("r =", round(cor(pausing_production_mat[, c(2,5)])[1,2], 2) )) +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 2]) + 0.1, label = "median") +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 2]) - 0.1, label = round(10^median(pausing_production_mat[, 2]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3, -2, -1, 0), labels=c(0.001, 0.01, 0.1, 1), limits = c(-3, 0)) +
  scale_y_continuous(name="Pausing Index", breaks=0:3, labels=10^(0:3), limits = c(-0.5, 3)) +
  labs(fill = "Genes", title = "2i")+
  theme(legend.position = "none")

g3 <- data.frame(pause_index = pausing_production_mat[, 3],
           RNA_production = pausing_production_mat[, 6]) %>%
  ggplot(aes(x = RNA_production, y = pause_index)) +
  geom_hex(binwidth=c(0.05, 0.08)) +
  geom_hline(yintercept = median(pausing_production_mat[, 3]), lty = 2, alpha = 0.5) +
  geom_text(x = -0.25, y = 3, label = paste("r =", round(cor(pausing_production_mat[, c(3,6)])[1,2], 2) )) +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 3]) + 0.1, label = "median") +
  geom_text(x = -0.2, y = median(pausing_production_mat[, 3]) - 0.1, label = round(10^median(pausing_production_mat[, 3]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3, -2, -1, 0), labels=c(0.001, 0.01, 0.1, 1), limits = c(-3, 0)) +
  scale_y_continuous(name="Pausing Index", breaks=0:3, labels=10^(0:3), limits = c(-0.5, 3)) +
  labs(title = "mTORi") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "FigS5_Pause_index_tx.png", path = "../figS5/figs",
       device = "png", width = 10, height = 3.5)

# plot gene body Pol II ~ RNA synthesis --------------------------------------------------
gene_body_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1.*Pol5p.*ALL.*bam$", full.names = T),
                           intervals = gene.gr, stranded = F, paired.end = "ignore")
gene_body_mat <- gene_body_mat / width(gene.gr) # convert to reads density

pol2s5p_input_size_factor <- readRDS("data/pol2s5p_input_size_factor.RData")
gene_body_mat <- sweep(gene_body_mat, 2, pol2s5p_input_size_factor, "/")

gene_body_mat <- gene_body_mat[, c(3, 1, 4)]
rownames(gene_body_mat) <- gene.gr$gene_id
colnames(gene_body_mat) <- c("Pol2_SL", "Pol2_2i", "Pol2_mTORi")
gene_body_mat <- log10(gene_body_mat)
gene_body_mat <- gene_body_mat[!apply(gene_body_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# add RNA synthesis
sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$gene_id %in% rownames(gene_body_mat), ]

sample_ord <- match(rownames(gene_body_mat), unique(sample_Tx_counts$gene_id))

gene_body_production_mat <- gene_body_mat
for ( i in c("SL", "2i_2d", "mTORi_2d")) {
  gene_body_production_mat <- cbind(gene_body_production_mat, # append nascent RNA production speed = copy * labeled rate
                                  rowSums(log10(sample_Tx_counts[sample_Tx_counts$Sample == i, c("Copy", "Labeled_rate")][sample_ord, ])) - log10(5))
}
gene_body_production_mat <- gene_body_production_mat[!apply(gene_body_production_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]
colnames(gene_body_production_mat)[4:6] <- c("tx_SL", "tx_2i", "tx_mTORi")
gene_body_production_mat[, 1:3] <- gene_body_production_mat[, 1:3] + 3 # arbitrary units of Pol II density

g4 <- ggplot(gene_body_production_mat, aes(x = tx_SL, y = Pol2_SL)) + # n = 5890
  geom_hex(bins= 50) +
  geom_hline(yintercept = median(gene_body_production_mat[, 1]), lty = 2, alpha = 0.5) +
  geom_text(x = 0.6, y = 2.8, label = paste("r =", round(cor(gene_body_production_mat[, c(1,4)])[1,2], 2) )) +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 1]) + 0.1, label = "median") +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 1]) - 0.1, label = round(10^median(gene_body_production_mat[, 1]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3:1), labels=c(0.001, 0.01, 0.1, 1, 10), limits = c(-3, 1)) +
  scale_y_continuous(name="Pol II S5p gene body (a.u.)", breaks=0:2, labels=10^(0:2), limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "SL")+
  theme(legend.position = "none")

g5 <- ggplot(gene_body_production_mat, aes(x = tx_2i, y = Pol2_2i)) + # n = 5890
  geom_hex(bins= 50) +
  geom_hline(yintercept = median(gene_body_production_mat[, 2]), lty = 2, alpha = 0.5) +
  geom_text(x = 0.6, y = 2.8, label = paste("r =", round(cor(gene_body_production_mat[, c(2,5)])[1,2], 2) )) +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 2]) + 0.1, label = "median") +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 2]) - 0.1, label = round(10^median(gene_body_production_mat[, 2]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3:1), labels=c(0.001, 0.01, 0.1, 1, 10), limits = c(-3, 1)) +
  scale_y_continuous(name="Pol II S5p gene body (a.u.)", breaks=0:2, labels=10^(0:2), limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "2i")+
  theme(legend.position = "none")

g6 <- ggplot(gene_body_production_mat, aes(x = tx_mTORi, y = Pol2_mTORi)) + # n = 5890
  geom_hex(bins= 50) +
  geom_hline(yintercept = median(gene_body_production_mat[, 3]), lty = 2, alpha = 0.5) +
  geom_text(x = 0.6, y = 2.8, label = paste("r =", round(cor(gene_body_production_mat[, c(3,6)])[1,2], 2) )) +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 3]) + 0.1, label = "median") +
  geom_text(x = 0.6, y = median(gene_body_production_mat[, 3]) - 0.1, label = round(10^median(gene_body_production_mat[, 3]), 1)) +
  theme_setting +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.95) +
  scale_x_continuous(name="\nNascent RNA (copy/min)", breaks=c(-3:1), labels=c(0.001, 0.01, 0.1, 1, 10), limits = c(-3, 1)) +
  scale_y_continuous(name="Pol II S5p gene body (a.u.)", breaks=0:2, labels=10^(0:2), limits = c(-0.5, 2.8)) +
  labs(fill = "Genes", title = "mTORi")+
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g4, g5, g6, ncol = 1),
       filename = "Fig4_Pol2S5p_genebody_tx.png", path = "figs",
       device = "png", width = 3.5, height = 10)

# plot RNA synthesis explained ~ features --------------------------------------------------
genes.intercect <- intersect.Vector(rownames(pausing_index_mat), rownames(gene_body_production_mat))

dat <- cbind(pausing_production_mat[genes.intercect, 1:3],
             gene_body_production_mat[genes.intercect, ],
             DHS = log(as.numeric(tss.gr$DHS_density[match(genes.intercect, tss.gr$gene_id)]))
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
g9 <- ggplot(R2_data, aes(x = Sample, y = R2, fill = Feature)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) +
  # coord_flip() +
  xlab("R-squared composition") +
  ylab("RNA synthesis explained") +
  scale_fill_manual(values = colors_9[c(2,8,1,5,4)]) +
  theme_setting +
  theme(legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))
ggsave(plot = g9, filename = "Fig4_RNA_synthesis_explained.png", path = "figs",
       device = "png", width = 4, height = 3)

# plot estimated speed ~ measured speed --------------------------------------------------
elongation_table <- readRDS("../fig5/data/elongation_table.RData")

dat_speed = data.frame(v_hat_SL = trim_quantile(dat$tx_SL - dat[, "Pol2_SL"]),
                       v_hat_2i = trim_quantile(dat$tx_2i - dat[, "Pol2_2i"]),
                       v_hat_mTORi = trim_quantile(dat$tx_mTORi - dat[, "Pol2_mTORi"])) 
dat_speed_cmp = data.frame(Sample = rep(c("SL", "2i", "mTORi"), each = nrow(dat_speed)),
                           Speed_hat = c(dat_speed$v_hat_SL, dat_speed$v_hat_2i, dat_speed$v_hat_mTORi))
rownames(dat_speed) <- rownames(dat)
dat_speed = dat_speed[!apply(dat_speed, 1, function(x) any(is.na(x))), ]
saveRDS(dat_speed, "data/samples_estimated_speed.RData")

gg_plot_speed <- function(x, y)
{
  ggplot(data.frame(x=x, y=y), aes(x = x, y = y)) +
    geom_hex(bins = 50) +
    scale_fill_viridis(option = "C", direction = -1, alpha = 0.8) +
    geom_text(x = -1.6, y = 1, label = paste("r =", round(cor(x, y), 3)), hjust = 0) +
    geom_text(x = -1.6, y = 0.85, label = paste("n =", length(x)), hjust = 0) +
    xlab("\nEstimated speed (a.u.)\n") +
    xlim(c(-3.5, -1)) +
    scale_y_continuous(name="Elongation speed (Kb/min)", breaks=(-1:1), labels=10^(-1:1), limits = c(-1.3, 1)) +
    theme_setting +
    theme(axis.ticks.x = element_blank(), 
          legend.position = "none") 
}
 
g7 <- gg_plot_speed(x = dat_speed$v_hat_SL, y = dat_speed$speed)

# plot estimated speed ~ samples ------------------------------------------------------------
library(ggridges)
dat_speed_cmp$Sample <- factor(dat_speed_cmp$Sample, levels = rev(c("SL", "2i", "mTORi")))
g8 <- ggplot(dat_speed_cmp, aes(x = Speed_hat, y = Sample, fill = Sample)) + 
  geom_density_ridges(rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2, lty = 2) +
  theme_setting +
  xlab("\nEstimated speed (a.u.)\n") +
  ylab("") +
  scale_fill_manual(values = colors_20[rev(c(13, 2, 7))]) +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none")
  

# plot RNA synthesis ~ estimated synthesis --------------------------------------------------

dat_mu <- data.frame(mu_hat = log10(elongation_table[genes.intercect, 2]) + dat[, "Pol2_SL"],
                     mu = trim_quantile(dat$tx_SL),
                     speed = log10(elongation_table[genes.intercect, 2]))
dat_mu <- dat_mu[!is.na(dat_mu$mu_hat), ]
g10 <- ggplot(dat_mu, aes(x = mu, y = mu_hat)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "A", direction = -1, alpha = 0.8) +
  geom_text(x = -0.2, y = 3, label = paste("r =", round(cor(dat_mu$mu_hat, dat_mu$mu), 3)), hjust = 0) +
  geom_text(x = -0.2, y = 2.6, label = paste("n =", nrow(dat_mu)), hjust = 0) +
  xlab("\nEstimated RNA synthesis (a.u.)\n") +
  ylab("RNA synthesis (copy/min)") +
  ylim(c(-1, 3)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none") 

# plot Pol II ~ speed -------------------------------------------------------
widths <- width(gene.gr[match(rownames(gene_body_mat), gene.gr$gene_id)])
dat_N_Pol2 <- data.frame(den = gene_body_mat,
                         Speed = log10(elongation_table[rownames(gene_body_mat), 2]))
dat_N_Pol2 = dat_N_Pol2[!apply(dat_N_Pol2, 1, function(x) any(is.na(x))), ]

g11 <- ggplot(dat_N_Pol2, aes(x = den.Pol2_SL, y = Speed)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "C", direction = -1, alpha = 0.8) +
  geom_text(x = -1.05, y = 1, label = paste("r =", round(cor(dat_N_Pol2$Speed, dat_N_Pol2$den.Pol2_SL), 3)), hjust = 0) +
  geom_text(x = -1.05, y = 0.8, label = paste("n =", nrow(dat_N_Pol2)), hjust = 0) +
  scale_y_continuous(name="Elongation speed (Kb/min)", breaks=(-1:1), labels=10^(-1:1), limits = c(-1.3, 1)) +
  xlab("\nPol II S5p density (a.u.)\n") +
  theme_setting +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none") 

# plot RNA synthesis ~ measured speed --------------------------------------------------

g12 <- ggplot(dat_mu, aes(x = mu, y = speed)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "A", direction = -1, alpha = 0.8) +
  geom_text(x = -0.15, y = 1, label = paste("r =", round(cor(dat_mu$speed, dat_mu$mu), 3)), hjust = 0) +
  geom_text(x = -0.15, y = 0.8, label = paste("n =", nrow(dat_mu)), hjust = 0) +
  scale_y_continuous(name="Elongation speed (Kb/min)", breaks=(-1:1), labels=10^(-1:1), limits = c(-1.3, 1)) +
  xlab("\nRNA synthesis (copy/min)\n") +
  theme_setting +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none") 

ggsave(plot = grid.arrange(g11, g12, g7, g10, g8, ncol = 2),
       filename = "Fig4_TTseq_Pol2S5p_speed_interplay.png", path = "figs",
       device = "png", width = 7, height = 10)

######################################################################################################
# Pol II elongation speed
######################################################################################################
source("F4_src_Tx_elongation.R")
# load processed data
pol2s5p_mat_sandwich_list <- readRDS("data/pol2s5p_gene_sandwich_mat.RData")
TT_gene_sandwich_mat_norm <- readRDS("data/TT_gene_sandwich_mat_norm.RData")
TT_Pol2s5p_sandwich_mat <- readRDS("data/TT_Pol2s5p_sandwich_mat.RData")

# plot average coverage
n_pos <- ncol(pol2s5p_mat_sandwich_list[[1]])

pol2s5p_data <- data.frame(position = rep(seq_len(n_pos), 3),
                           value = c(colMeans(log1p(pol2s5p_mat_sandwich_list[[1]])),
                                     colMeans(log1p(pol2s5p_mat_sandwich_list[[2]])),
                                     colMeans(log1p(pol2s5p_mat_sandwich_list[[3]]))),
                           sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                           levels = c("SL", "2i", "mTORi"))
                           )

ttseq_data <- data.frame(position = rep(seq_len(n_pos), 3),
                           value = c(colMeans(log1p(TT_gene_sandwich_mat_norm[[1]])),
                                     colMeans(log1p(TT_gene_sandwich_mat_norm[[2]])),
                                     colMeans(log1p(TT_gene_sandwich_mat_norm[[3]]))),
                           sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                     levels = c("SL", "2i", "mTORi"))
)

TT_Pol2s5p_data <- data.frame(position = rep(seq_len(n_pos), 3),
                         value = c(colMeans(TT_Pol2s5p_sandwich_mat[[1]]),
                                   colMeans(TT_Pol2s5p_sandwich_mat[[2]]),
                                   colMeans(TT_Pol2s5p_sandwich_mat[[3]])),
                         sample = factor(c(rep("SL", n_pos), rep("2i", n_pos), rep("mTORi", n_pos)),
                                         levels = c("SL", "2i", "mTORi"))
)

g1.1 <- ggplot(pol2s5p_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 220, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 2) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Pol II-S5p)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(20, 220), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "top")

g1.2 <- ggplot(ttseq_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 220, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 2) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Nascent RNA)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(20, 220), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "none")

g1.3 <- ggplot(TT_Pol2s5p_data, aes(x = position, y = value, color = sample)) +
  geom_rect(aes(xmin = -Inf, xmax = 20, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 220, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(size = 2) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  xlab("") +
  ylab("log(Nascent RNA / Pol II-S5p)") +
  labs(color = "Sample") +
  scale_x_continuous(breaks = c(20, 220), labels = c("TSS", "TTS")) +
  theme_setting +
  theme(axis.text.x = element_text(size=14),
        legend.position = "none")

ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, ncol = 1, heights = c(4.7, 4, 4)),
       filename = "Fig4_TTseq_Pol2S5p_coverage.png", path = "figs",
       device = "png", width = 3.5, height = 10)
