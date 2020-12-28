#-----------------------------------------------------------------------------------------------------
# H3K4me3 and H3K27me3 changes correlation test with tx, half-life and copy number changes
#-----------------------------------------------------------------------------------------------------

# MINUTE ChIP Histone mark changes --------------------------------------------------------------------------------------------
if (T) {
  bw_files <- list.files('/mnt/E0767589767560E8/UPPMAX/SL2i', pattern = 'H3|H4.*bw', full.names = T)
  
  design_mat <- data.frame(condition = gsub("(.*_)(.*)_(CTR|INH).*", "\\2", bw_files),
                           mark = gsub(".*SL2i/(H.*)_(SL|2i).*", "\\1", bw_files),
                           treat = gsub("(.*_)(.*)_.*", "\\2", bw_files))
  # get ChIP RPK
  histone_counts <- .countBW(bw_files = bw_files,
                             intervals = TU.DE.mm9.gr + 500,
                             fast = T) %>%
                    trim_quantile() %>%
                    as.data.frame()
  
  colnames(histone_counts) <- do.call(sprintf,  c(design_mat[c(2,1,3)], '%s_%s_%s'))
  
  histone_tss_counts <- .countBW(bw_files = bw_files,
                                 intervals = TU.TSS.mm9.gr,
                                 fast = T) %>%
                        trim_quantile() %>%
                        as.data.frame()

  colnames(histone_tss_counts) <- do.call(sprintf, c(design_mat[c(2,1,3)], '%s_%s_%s'))
  
  fold_changes <- function(X, Y)
    ifelse( ((X*Y < 1e-5) | abs(X - Y) < 1e-18 ), NA, 1) * (X - Y)
    
  rlog_transform <- function(dat, .sample) {
    dds <- DESeqDataSetFromMatrix(round(dat), 
                                  colData = data.frame(Sample = .sample), 
                                  design = ~Sample)
    sizeFactors(dds) <- 1
    dds_histone %>% rlog() %>% assay() %>% as.data.frame() %>% return()
  }
  
  # stablize low count genes with rlog transformation to balance log2FC distribution
  require(DESeq2)
  histone_counts <- rlog_transform(histone_counts, design_mat$mark)
  histone_tss_counts <- rlog_transform(histone_tss_counts, design_mat$mark)
  
  histone_counts[, grepl("H3K4m3", bw_files)] <- histone_tss_counts[, grepl("H3K4m3", bw_files)]
  
  # make the mean and change table
  mark_cmp_2i <- c("H3K4m3", "H3K27m2", "H3K27m3")
  mark_cmp_mTORi <- c("H3K4m3", "H3K27m2", "H3K27m3", "H3K27ac", "H4K5ac", "H4K12ac")
  
  res_mat <- cbind("H3K4m3" = rowMeans(histone_counts[, design_mat$mark == "H3K4m3"]),
                   "H3K27m2" = rowMeans(histone_counts[, design_mat$mark == "H3K27m2"]),
                   "H3K27m3" = rowMeans(histone_counts[, design_mat$mark == "H3K27m3"]),
                   "H3K27ac" = rowMeans(histone_counts[, design_mat$mark == "H3K27ac"]),
                   "H4K5ac" = rowMeans(histone_counts[, design_mat$mark == "H4K5ac"]),
                   "H4K12ac" = rowMeans(histone_counts[, design_mat$mark == "H4K12ac"]),
                   
                   fold_changes(histone_counts[, with(design_mat, 
                                                       treat == "CTR" &
                                                         condition == "2i" &
                                                         mark %in% mark_cmp_2i)],
                                 histone_counts[, with(design_mat, 
                                                       treat == "CTR" &
                                                         condition == "SL" &
                                                         mark %in% mark_cmp_2i)]),
                   
                   fold_changes(histone_counts[, with(design_mat,
                                                      treat == "INH" &
                                                        condition == "SL" &
                                                        mark %in% mark_cmp_mTORi)],
                                 histone_counts[, with(design_mat,
                                                       treat == "CTR" &
                                                         condition == "SL" &
                                                         mark %in% mark_cmp_mTORi)]))
}

# link epi-mark changes to transcription changes
dat_h_tx_fc <- data.frame("H3K4m3_2i_FC" = res_mat$H3K4m3_2i_CTR,
                          "H3K27m3_2i_FC" = res_mat$H3K27m3_2i_CTR,
                          "H3K4m3_mTORi_FC" = res_mat$H3K4m3_SL_INH,
                          "H3K27m3_mTORi_FC" = res_mat$H3K27m3_SL_INH,
                          "LRNA_FC_2i" = TU.DE.mm9.gr$log2FoldChange_LRNA_sp_2i,
                          "LRNA_FC_mTORi" = TU.DE.mm9.gr$log2FoldChange_LRNA_sp_mTORi)

plot_FC <- function(dat, .xlab, .ylab, .col, .title, .xlim = c(-5, 5)) {
  dat <- dat %>%
    `colnames<-`(c("x", "y")) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate("Density" = get_dens(x, y, n.grid = 100))
  
  ggplot(dat, aes(x = x, y = y, color = Density)) +
    geom_point(pch = 19, size = 0.3) +
    # geom_hex(aes(color = Density), bins = 60) +
    geom_text(x = .xlim[2] - (.xlim[2] - .xlim[1])/3, y = 4.7, 
              label = paste("r =", round(cor(dat$x, dat$y), 3)), 
              hjust = 0, vjust = 1, color = "black") +
    scale_color_viridis_c(option = .col, direction = -1, end = 0.8) +
    xlim(.xlim) + ylim(c(-5, 5)) +
    xlab(.xlab) + ylab(.ylab) + ggtitle(.title) +
    theme_setting +
    theme(legend.position = "none") 
}

# 2i
g_list <- list()
for (i in seq_along(unique(umap_dat$Cluster)) ) {
  g_list <- c(g_list, 
              list(plot_FC(dat_h_tx_fc[tmp.idx[umap_dat$Cluster == i], 
                                       c("H3K4m3_2i_FC", "LRNA_FC_2i")], 
                           "H3K4m3 2i Log2FC", ifelse(i == 1, "Tx 2i Log2FC", ""), "C",
                           paste0("C", i), c(-3, 2))))
}

for (i in seq_along(unique(umap_dat$Cluster)) ) {
  g_list <- c(g_list,
              list(plot_FC(dat_h_tx_fc[tmp.idx[umap_dat$Cluster == i], 
                                       c("H3K27m3_2i_FC", "LRNA_FC_2i")], 
                           "H3K27m3 2i Log2FC ", ifelse(i == 1, "Tx 2i Log2FC", ""), "D", 
                           "", c(-3, 4))))
}
ggsave(plot = do.call(grid.arrange, c(g_list, nrow = 2)),
       filename = paste0("Fig3_2i_TU_tx_log2FC_H3K4me3_H3K27me3_log2FC.png"), 
       path = "../fig3/figs/",
       device = "png", width = 18, height = 6)

# mTORi
g_list <- list()
for (i in seq_along(unique(umap_dat$Cluster)) ) {
  g_list <- c(g_list, 
              list(plot_FC(dat_h_tx_fc[tmp.idx[umap_dat$Cluster == i], 
                                       c("H3K4m3_mTORi_FC", "LRNA_FC_mTORi")], 
                           "H3K4m3 mTORi Log2FC", ifelse(i == 1, "Tx mTORi Log2FC", ""), "C",
                           paste0("C", i), c(-3, 2))))
}
for (i in seq_along(unique(umap_dat$Cluster)) ) {
  g_list <- c(g_list,
              list(plot_FC(dat_h_tx_fc[tmp.idx[umap_dat$Cluster == i], 
                                       c("H3K27m3_mTORi_FC", "LRNA_FC_mTORi")], 
                           "H3K27m3 mTORi Log2FC ", ifelse(i == 1, "Tx mTORi Log2FC", ""), "D", 
                           "", c(-3, 4))))
}
ggsave(plot = do.call(grid.arrange, c(g_list, nrow = 2)),
       filename = paste0("Fig3_mTORi_TU_tx_log2FC_H3K4me3_H3K27me3_log2FC.png"), 
       path = "../fig3/figs/",
       device = "png", width = 18, height = 6)


