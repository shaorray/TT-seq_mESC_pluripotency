

# plot comparison  -------------------------------------------------------------------------
simple_cor_plot <- function(x, y, viridis_type = "A", color_end = 1, .xlab, .ylab)
{
  data.frame(x = log(x), y = log(y)) %>%
    dplyr::filter(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)) %>%
    dplyr::filter(x > quantile(x, 0.05) & y > quantile(y, 0.05)) %>%
    dplyr::filter(x < quantile(x, 0.99) & y < quantile(y, 0.99)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_hex(bins= 50) +
    scale_fill_viridis_c(option = viridis_type, direction = -1, end = color_end) +
    xlab(.xlab) + ylab(.ylab) +
    theme_setting +
    theme(legend.position = "none")
}

# load data  -------------------------------------------------------------------------------
sample_conditions <- c("SL", "2i_2d", "mTORi_1d")
L_RNA_sf <- readRDS("../data/LRNA.sizefactor.RData")
gene.gr <- readRDS("data/gene.gr.RData") # 
gene.gr <- gene.gr[width(gene.gr) > 2000 & width(gene.gr) < 1000000]

# ChIP coverage ----------------------------------------------------------------------------
# estimate MINUTE ChIP sample sizes
pol2s5p_input_size_factor = .countBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                                             pattern = "P1_IN_.*ALL.*bam$", full.names = T),
                                      gene.gr) %>% SizeFactorCal
saveRDS(pol2s5p_input_size_factor, "data/pol2s5p_input_size_factor.RData")

# read Pol2s5p coverage
bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                       pattern = "P1.*Pol.*ALL.*bam$", full.names = T)

pol2s5p_mat_sandwich <- readBam(bam_files = bam_files,
                                intervals = gene.gr,
                                pair_end = T,
                                stranded = F,
                                flanks = c(2000, 4000),
                                new_lens = c(20, 200, 40))

# normalise ChIP coverage between samples
pol2s5p_mat_sandwich_list <- list()
for (i in c(3, 1, 4)) {
  pol2s5p_mat_sandwich_list <- c(pol2s5p_mat_sandwich_list,
                                 list(pol2s5p_mat_sandwich[[i]] / pol2s5p_input_size_factor[i]))
}
names(pol2s5p_mat_sandwich_list) <- sample_conditions
saveRDS(pol2s5p_mat_sandwich_list, "data/pol2s5p_gene_sandwich_mat.RData")

# TTseq coverage  ----------------------------------------------------------------------------
bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                       pattern = "LRNA.*(SL_|_2i_2d|mTORi_).*bam$", full.names = T)

TT_gene_sandwich_mat <- readBam(bam_files = bam_files,
                                intervals = gene.gr,
                                pair_end = T,
                                stranded = F, # include upstream reads on the main gene coverage
                                flanks = c(2000, 4000),
                                new_lens = c(20, 200, 40))

names(TT_gene_sandwich_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1", names(TT_gene_sandwich_mat))

# normalise TTseq coverage by spike-in size factors
for (i in names(TT_gene_sandwich_mat)) {
  TT_gene_sandwich_mat[[i]] <- TT_gene_sandwich_mat[[i]] / L_RNA_sf[i]
}
# averaging replicates
TT_gene_sandwich_mat_norm <- list()
for (i in c("SL", "2i", "mTORi")) {
  idx <- grep(i, names(TT_gene_sandwich_mat))
  TT_gene_sandwich_mat_norm <- c(TT_gene_sandwich_mat_norm,
                                 list(Reduce("+", TT_gene_sandwich_mat[idx]) / length(idx)))
}

names(TT_gene_sandwich_mat_norm) <- sample_conditions
saveRDS(TT_gene_sandwich_mat_norm, "data/TT_gene_sandwich_mat_norm.RData")

TT_Pol2s5p_sandwich_mat <- list()
for (i in sample_conditions) {
  TT_Pol2s5p_sandwich_mat <- c(TT_Pol2s5p_sandwich_mat,
                          list(log1p(TT_gene_sandwich_mat_norm[[i]]) - log1p(pol2s5p_mat_sandwich_list[[i]])) )
}
names(TT_Pol2s5p_sandwich_mat) <- sample_conditions

saveRDS(TT_Pol2s5p_sandwich_mat, "data/TT_Pol2s5p_sandwich_mat.RData")

