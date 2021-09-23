#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# read in kallisto tx counts on gencode.vM17.annotation and spike-in RNAs
# and save reads counts for normalization
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

# load kallisto counts, mm10
filenames <- sort(list.files('../data/kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                     as = 'matrix', what = "tpm")
colnames(count_table) <- sampleNewName

# split to GENCODE transcripts, annotated TUs, spike-in RNAs
txRPK <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T, is_gene_sum = F)

tuRPK <- count_table[!grepl("^chrS|^ENS", rownames(count_table)), ]
tuRPK <- tuRPK[rowSums(tuRPK) > 0, ]

spRPK <- count_table[grep("^chrS", rownames(count_table)), ]

saveRDS(txRPK, "../data/txRPK_SL_2i.RData") # raw tmp without normalization
saveRDS(tuRPK, "../data/tuRPK_SL_2i.RData") # mm10 TU annotation kallisto tpm
saveRDS(spRPK, "../data/spRPK_SL_2i.RData")

# matrix for spike-in normalization
sp_F_mat <- spRPK[, grepl("FRNA", colnames(spRPK))]
colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
FRNA.sizefactor <- SizeFactorCal(sp_F_mat)

sp_L_mat <- spRPK[, grepl("LRNA",colnames(spRPK))]
colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins

saveRDS(FRNA.sizefactor, "../data/FRNA.sizefactor.RData")
saveRDS(LRNA.sizefactor, "../data/LRNA.sizefactor.RData")

if (F) { # spike-in size factors for read count normalization
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table) <- sampleNewName
  
  # norm to RPK
  eff_lengths <- SummarizedExperiment::readKallisto(paste0(filenames[1], "/abundance.tsv"),
                                                    as = 'matrix', what = "eff_length")
  count_table <- count_table / c(eff_lengths) * 1e3
  
  spRPK <- count_table[grep("^chrS", rownames(count_table)), ]
  
  sp_F_mat <- spRPK[, grepl("FRNA", colnames(spRPK))]
  colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
  FRNA.sizefactor <- SizeFactorCal(sp_F_mat)
  
  sp_L_mat <- spRPK[, grepl("LRNA",colnames(spRPK))]
  colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
  LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins
  
  saveRDS(FRNA.sizefactor, "../data/FRNA.sizefactor.RC.RData")
  saveRDS(LRNA.sizefactor, "../data/LRNA.sizefactor.RC.RData")
}

## spike-ins table, for labelel rate estimation
spikein_lens <- c("chrS2" = 1.023, "chrS4" = 1.033, "chrS5" = 1.042,
                  "chrS8" = 1.124, "chrS9" = 1.061, "chrS12" = 1.124) 
# convert RPK to abundance by multiplying spikein lengths, since spikeins were mixed by weight
SampleSpCounts <- NULL
for(i in seq_len(ncol(sp_F_mat))){
  SampleSpCounts <- rbind(SampleSpCounts,
                          data.frame(FRNA = sp_F_mat[, i] / FRNA.sizefactor[i] * spikein_lens,
                                     LRNA = sp_L_mat[, i] / LRNA.sizefactor[i] * spikein_lens,
                                     Sample = colnames(sp_L_mat)[i],
                                     SpikeIns = rownames(sp_F_mat),
                                     W = c(1, 0.1, 1, 0.1, 1, 0.1),
                                     R = c(1, 1, 0.1, 0.1, 0, 0) ))
}
saveRDS(SampleSpCounts, "../data/SampleSpikeCounts.RData")

# normalise tx tpm with spike-in RNA size factor, for parameter estimation
tx_F_mat <- txRPK[, grep("FRNA", colnames(txRPK))]
colnames(tx_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tx_F_mat))
tx_F_mat <- sweep(tx_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tx_L_mat <- txRPK[, grep("LRNA", colnames(txRPK))]
colnames(tx_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tx_L_mat))
tx_L_mat <- sweep(tx_L_mat, 2, LRNA.sizefactor, '/')

saveRDS(tx_L_mat, "../data/tx_L_mat.RData")
saveRDS(tx_F_mat, "../data/tx_F_mat.RData")

# normalise non-coding TU tpm with spike-in RNA size factor
tu_F_mat <- tuRPK[, grep("FRNA", colnames(tuRPK))]
colnames(tu_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tu_F_mat))
tu_F_mat <- sweep(tu_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tu_L_mat <- tuRPK[, grep("LRNA", colnames(tuRPK))]
colnames(tu_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tu_L_mat))
tu_L_mat <- sweep(tu_L_mat, 2, LRNA.sizefactor, '/')

saveRDS(tu_F_mat, "../fig1/data/tu_FRNA_RPK_norm.RData")
saveRDS(tu_L_mat, "../fig1/data/tu_LRNA_RPK_norm.RData")

# process mm9 tx counts ------------------------------------------------------------------------------------------
if (F) {
  # count with combined reference GENCODE vM20 and ncRNA annotation
  filenames <- sort(list.files('../data/kallisto_output_mm9/', full.names = T)) 
  sampleNewName <- gsub(".*/", "\\2", filenames)
  
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table) <- sampleNewName
  
  txRC <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
    keepOneTx(rowname_gene_id = T, is_gene_sum = T)
  txRC <- txRC[rowSums(txRC) > 0, ]
  
  # save read count
  saveRDS(txRC, "../data/txRC_SL_2i_mm9.RData") # raw count without normalization

  # matrix for spike-in normalization
  sp_F_mat <- spRC[, grepl("FRNA", colnames(spRC))]
  colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
  FRNA.sizefactor <- SizeFactorCal(sp_F_mat)
  
  sp_L_mat <- spRC[, grepl("LRNA",colnames(spRC))]
  colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
  LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins
  
  saveRDS(list("FRNA.sizefactor" = FRNA.sizefactor, 
               "LRNA.sizefactor" = LRNA.sizefactor),
          "../fig1/data/sizefactor.list.mm9.RData")
  
  # normalise with spike-ins size factors
  tx_F_mat <- txRC[, grepl("FRNA", colnames(txRC))]
  colnames(tx_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tx_F_mat))
  tx_F_mat <- sweep(tx_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors
  
  tx_L_mat <- txRC[, grepl("LRNA", colnames(txRC))]
  colnames(tx_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tx_L_mat))
  tx_L_mat <- sweep(tx_L_mat, 2, LRNA.sizefactor, '/')
  
  # saveRDS(list("tx_F_mat" = tx_F_mat, "tx_L_mat" = tx_L_mat, 
  #              "tu_F_mat" = tu_F_mat, "tu_L_mat" = tu_L_mat), 
  #         "../fig1/data/tx_tu_RPK_norm.mm9.RData")
  
  library(DESeq2)
  dds_FRNA <- DESeqDataSetFromMatrix(round(as.matrix(tx_F_mat)),
                                     colData = data.frame(condition = gsub("(.*)_rep.", "\\1", colnames(tx_F_mat)) ),
                                     design = ~ condition)
  dds_LRNA <- DESeqDataSetFromMatrix(round(as.matrix(tx_L_mat)),
                                     colData = data.frame(condition = gsub("(.*)_rep.", "\\1", colnames(tx_F_mat)) ),
                                     design = ~ condition)
  dds_FRNA <- DESeq(dds_FRNA)
  dds_LRNA <- DESeq(dds_LRNA)
  
  res_FRNA_2i <- results(dds_FRNA, contrast = c("condition", "2i_2d", "SL"))
  res_FRNA_mTORi <- results(dds_FRNA, contrast = c("condition", "mTORi_1d", "SL"))
  res_FRNA_SL2i <- results(dds_FRNA, contrast = c("condition", "SL2i_2d", "SL"))
  
  res_LRNA_2i <- results(dds_LRNA, contrast = c("condition", "2i_2d", "SL"))
  res_LRNA_mTORi <- results(dds_LRNA, contrast = c("condition", "mTORi_1d", "SL"))
  res_LRNA_SL2i <- results(dds_LRNA, contrast = c("condition", "SL2i_2d", "SL"))
  
  res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                         keys = rownames(res_FRNA_2i),
                         keytype = "GENEID",
                         columns = "GENENAME")
  
  # Bulut et al. 2016 2i marker genes are specific to their own measurements 
  # "Spi1", "Prdm16", "Bmp7", "Sp100", "Dazl", "Trpm1", "Crxos"
  marker_gene_2i <- c("Myc", "Zic3", "Nanog", "Utf1", "Dnmt3l", "Etv4", "Id1", "Lefty1", 
                      "Tfcp2l1", "Fgf10", "Cdh2", "Lefty2", "Zic1", "Neurog2", "Sox1", "Sox17")
  marker_gene_id_2i <- res$GENEID[match(marker_gene_2i, res$SYMBOL)]
  
  marker_gene_mTORi <- c("Nphs1", "Hbp1", "Kirrel2", "Platr7", "Lefty1", "Txnip",
                         "Pdcd4", "Myrf", "Zfp652", "Aplp1", "Pim3", "Meg3", "Lefty2", "Grhl2")
  
  dat_log2FC <- data.frame(log2FC_2i = c(res_FRNA_2i[marker_gene_id_2i, "log2FoldChange"],
                                         res_LRNA_2i[marker_gene_id_2i, "log2FoldChange"]),
                           
                           log2FC_SL2i = c(res_FRNA_SL2i[marker_gene_id_2i, "log2FoldChange"],
                                           res_LRNA_SL2i[marker_gene_id_2i, "log2FoldChange"]),
                           Type = c(rep("FRNA", each = length(marker_gene_2i)),
                                    rep("LRNA", each = length(marker_gene_2i))), 
                           gene_name = rep(marker_gene_2i, 2),
                           gene_label = c(marker_gene_2i, rep("", length(marker_gene_2i))))
  
  ggplot(dat_log2FC, aes(x = log2FC_2i, y = log2FC_SL2i,
                         shape = Type, label = gene_label,
                         color = gene_name, group = gene_name)) +
    geom_abline(intercept = 0, slope = 1, color = add.alpha("grey50", 0.5)) +
    geom_vline(xintercept = 0, color = add.alpha("grey50", 0.5)) +
    geom_hline(yintercept = 0, color = add.alpha("grey50", 0.5)) +
    geom_line(color = "grey50", lty = 2, size = 0.5) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(1, 4)) +
    ggtitle("DESeq2 of 2i marker genes (n=16)") +
    xlab("log2FC 2i 2d") + ylab("log2FC SL2i 2d") +
    geom_text(size = 4, hjust = -0.2, vjust = 0.7, check_overlap = T) +
    guides(color = FALSE) +
    theme_setting
    
  ggsave(filename = "FigS2_dot_plot_2i_marker_gene_2i_vs_SL2i.png",
         path = "../figS2/figs", width = 5, height = 4, device = "png")
  
  
  g1 <- plot_scatter(dat = data.frame(x = res_FRNA_2i[names(gene.gr), "log2FoldChange"],
                                      y = res_FRNA_SL2i[names(gene.gr), "log2FoldChange"]), 
                     .xlab = "2i 2d", .ylab = "SL2i 2d", 
                     xlim = c(-5, 5), ylim = c(-5, 5)) + ggtitle("FRNA log2FoldChange")
  
  g2 <- plot_scatter(dat = data.frame(x = res_LRNA_2i[names(gene.gr), "log2FoldChange"],
                                      y = res_LRNA_SL2i[names(gene.gr), "log2FoldChange"]), 
                     .xlab = "2i 2d", .ylab = "SL2i 2d", 
                     xlim = c(-5, 5), ylim = c(-5, 5)) + ggtitle("LRNA log2FoldChange")
  
  ggsave(plot = grid.arrange(g1, g2, nrow = 1),
         filename = "FigS2_dot_plot_gene_log2FC_2i_vs_SL2i.png",
         path = "../figS2/figs", width = 6, height = 3, device = "png")
}

if (F) { # external mESC TT-seq data
  TT_gene_RPK <- .countBam(bam_files = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA_mm9/2018_Jan_TTseq_TX1072_1/bam",
                                                  "*bam$", full.names = T),
                           intervals = gene.gr, 
                             stranded = T, paired.end = "ignore") / width(gene.gr) * 1e3
  
  WT26_gene_RPK <- .countBam(bam_files = "/mnt/0E471D453D8EE463/TT_seq_data/161226_H33KO_ATRX_DAXX_KO/bam_mm9/TTSEQ-WT26.Aligned.sortedByCoord.out.bam",
                             intervals = gene.gr, 
                             stranded = T, paired.end = "ignore") / width(gene.gr) * 1e3
  
  H33WT_gene_RPK <- .countBam(bam_files = "/mnt/0E471D453D8EE463/TT_seq_data/161226_H33KO_ATRX_DAXX_KO/bam_mm9/TTSEQ-H33WT.Aligned.sortedByCoord.out.bam",
                              intervals = gene.gr, 
                              stranded = T, paired.end = "ignore") / width(gene.gr) * 1e3
  
  TT_E14_SL_RPK <- .countBam(bam_files = "/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9/LRNA_SL_rep1.Aligned.sortedByCoord.out.bam",
                              intervals = gene.gr, 
                              stranded = T, paired.end = "ignore") / width(gene.gr) * 1e3
  
  dat_TT <- data.frame(Jan_TTseq = TT_gene_RPK[, 1], 
                       E14_SL = TT_E14_SL_RPK[, 1],
                       WT26 = WT26_gene_RPK[, 1],
                       H33WT = H33WT_gene_RPK[, 1]) %>% 
    log2() %>%
    as.data.frame() %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & Jan_TTseq > (-5) & E14_SL > (-5))
  
  g1 <- ggplot((dat_TT), aes(x = E14_SL, y = Jan_TTseq, 
                     color = get_dens(E14_SL, Jan_TTseq))) +
    geom_point(cex = 0.2) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2, 
             label = paste0("     r = ", round(cor(dat_TT, method = "spearman")[2, 1], 3), 
                            "\nn = ", nrow(dat_TT))) +
    xlab("log2 Mean RPK E14 SL") + ylab("log2 RPK Jan et al. 0h") + 
    ggtitle("TT-seq WT mESC SL") +
    scale_color_viridis_c(option = "C", direction = -1, begin = 0.05, end = 0.9) +
    theme_setting +
    theme(legend.position = "none")
  
  g2 <- ggplot(dat_TT, aes(x = E14_SL, y = WT26, 
                       color = get_dens(E14_SL, WT26))) +
    geom_point(cex = 0.2) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2, 
             label = paste0("     r = ", round(cor(dat_TT, method = "spearman")[2, 3], 3), 
                            "\nn = ", nrow(dat_TT))) +
    xlab("log2 Mean RPK E14 SL") + ylab("log2 RPK WT26") + 
    ggtitle("TT-seq WT mESC SL") +
    scale_color_viridis_c(option = "C", direction = -1, begin = 0.05, end = 0.9) +
    theme_setting +
    theme(legend.position = "none")
  
  g3 <- ggplot(dat_TT, aes(x = E14_SL, y = H33WT,
                       color = get_dens(E14_SL, H33WT))) +
    geom_point(cex = 0.2) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2, 
             label = paste0("     r = ", round(cor(dat_TT, method = "spearman")[2, 4], 3), 
                            "\nn = ", nrow(dat_TT))) +
    xlab("log2 Mean RPK E14 SL") + ylab("log2 RPK H33WT") + 
    ggtitle("TT-seq WT mESC SL") +
    scale_color_viridis_c(option = "C", direction = -1, begin = 0.05, end = 0.9) +
    theme_setting +
    theme(legend.position = "none")
  
  ggsave(grid.arrange(g1, g2, g3, nrow = 1), 
         filename = "../figS1/figs/FigS1_Scatter_mESC_TTseq_comparison.png",
         width = 12, height = 4)
}
# ------------------------------------------------------------------------------------------
meanSampleCounts <- function(tx_mat)
{ # this function averaging replicates
  tmp_mat = NULL
  for( i in unique(colnames(tx_mat)) )
  {
    idx = colnames(tx_mat) == i

    if( sum(idx) > 1)
    {
      tmp_mat <- cbind(tmp_mat, rowMeans(tx_mat[, idx]))
    } else {
      tmp_mat <- cbind(tmp_mat, tx_mat[, idx])
    }
  }
  colnames(tmp_mat) <- unique(colnames(tx_mat))
  return(tmp_mat)
}

# ---------------------------------------------------------------------------- #
# PCA plots of public data
if (F) {
  # load kallisto counts, SL vs 2i public data
  filenames <- sort(list.files('../data/kallisto_output_SL_2i/', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
  sampleNewName <- gsub(".*/", "\\2", filenames) %>% 
    gsub("201._", "", .) %>% 
    gsub("RNAseq_|RNASeq_|RNA-Seq_|_RNA-seq|_RNA-Seq", "", .) %>%
    gsub("mES_WT_|E14.|ESC", "", .)

  count_table_SL2i <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table_SL2i) <- sampleNewName
  
  txRC_SL2i <- count_table_SL2i[grepl("^ENS", rownames(count_table_SL2i)), ] %>%
    keepOneTx(rowname_gene_id = T, is_gene_sum = T)
  
  txRC_SL2i <- sapply(unique(gsub("_.$|_rep.", "", sampleNewName)), 
                      function(x) {
                        tmp <- txRC_SL2i[, grep(x, colnames(txRC_SL2i))]
                        if (!is.null(dim(tmp))) {
                          rowMeans(tmp)
                        } else {
                          tmp
                        }
                      })
  txRC_SL2i_LFC <- data.frame("Galonska_2i_24h" = txRC_SL2i[, 1] / txRC_SL2i[, 3],
                          "Galonska_2i_3d" = txRC_SL2i[, 2] / txRC_SL2i[, 3],
                          "Bulut_2i" = txRC_SL2i[, 4] / txRC_SL2i[, 6],
                          "Bulut_mTORi" = txRC_SL2i[, 5] / txRC_SL2i[, 6],
                          "Bulut_v6.5_2i" = txRC_SL2i[, 7] / txRC_SL2i[, 9],
                          "Bulut_v6.5_mTORi" = txRC_SL2i[, 8] / txRC_SL2i[, 9],
                          "Finley_2i" = txRC_SL2i[, 10] / txRC_SL2i[, 11],
                          "Joshi_2i" = txRC_SL2i[, 12] / txRC_SL2i[, 13],
                          "Marks_2i" = txRC_SL2i[, 14] / txRC_SL2i[, 15]) %>% log2()
  
  txRC_FRNA <- readRDS("../data/txRC_SL_2i_mm9.RData") 
  txRC_FRNA <- txRC_FRNA[, grepl("LRNA", colnames(txRC_FRNA))]
  txRC_FRNA <- sapply(unique(gsub("_rep.", "", colnames(txRC_FRNA))), 
                      function(x) {
                        tmp <- txRC_FRNA[, grep(x, colnames(txRC_FRNA))]
                        if (!is.null(dim(tmp))) {
                          rowMeans(tmp)
                        } else {
                          tmp
                        }
                      })
  txRC_FRNA_LFC <- data.frame("2i_2d" = txRC_FRNA[,1] / txRC_FRNA[,5], 
                              # "2i_7d" = txRC_FRNA[,2] / txRC_FRNA[,5], # no replicate
                              "SL2i_2d" = txRC_FRNA[,6] / txRC_FRNA[,5], 
                              "mTORi_1d" = txRC_FRNA[,3] / txRC_FRNA[,5], 
                              "mTORi_2d" = txRC_FRNA[,4] / txRC_FRNA[,5])
  gene.ov <- intersect.Vector(rownames(txRC_FRNA_LFC), rownames(txRC_SL2i_LFC)) 
  
  txRC_all_LFC <- cbind(txRC_FRNA_LFC[gene.ov, ],
                    txRC_SL2i_LFC[gene.ov, ])
  txRC_all_LFC <- txRC_all_LFC[is.finite(rowSums(txRC_all_LFC)) & !is.na(rowSums(txRC_all_LFC)), ]
 
  pca_all_LFC <- prcomp(t(txRC_all_LFC[intersect.Vector(res$GENEID[res$SYMBOL %in% unlist(pluripotent_markers)], # a list from "FigS2_RNA_turnover_comparison.R"
                                                        rownames(txRC_all_LFC)), ]))
  
  plot_pca <- function(pca, sample_names) {
    pca$sdev <- pca$sdev^2
    pc_var <- pca$sdev[1:2] / sum(pca$sdev) * 100
    pc <- pca$x[, 1:2] %>% as.data.frame()
    pc$sample_names <- sample_names
    pc$Sample <- gsub("_.*", "", sample_names)
    pc$Group <- c("2i", "mTORi")[grepl("2i", sample_names) + 1]
    
    hull <- pc %>%
      group_by(Group) %>%
      slice(chull(PC1, PC2))
    
    ggplot(pc, aes(x = PC1, y = PC2, color = Group, group = Group)) +
      geom_point(size = 4) +
      ggrepel::geom_text_repel(aes(label = sample_names), size = 4) +
      geom_polygon(data = hull, alpha = 0 ) +
      # ggforce::geom_mark_ellipse() +
      # stat_ellipse(aes(x = PC1, y = PC2, group = Group), type = "norm") +
      # scale_color_manual(values = colors_20[c(2, 16, 20, 7, 13, 10)]) +
      xlim(range(pc$PC1) * 1.5 ) + 
      xlab(paste0("PC1 ", round(pc_var[1]), "% variance")) +
      ylab(paste0("PC2 ", round(pc_var[2]), "% variance")) +
      theme_setting + 
      theme(legend.position = "none")
  }
  
  plot_pca(pca = pca_all_LFC, sample_names = gsub("^X", "", colnames(txRC_all_LFC))) + 
    ggtitle("Pluripotent genes log2FC (n=36)")
  
  ggsave(filename = "FigS1_PCA_2i_mTORi_log2FC_public_data_comparison.png", 
         path = "../figS1/figs",
         device = "png", width = 5, height = 5)
}
