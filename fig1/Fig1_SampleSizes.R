# measure nascent transcription and steady state RNA amount on cell level
# Rui Shao
# 2019 Sep

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("F1_src_LoadReadCounts.R")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

tx_FRNA_sf <- tx_F_mat %>% SizeFactorCal()
tx_LRNA_sf <- tx_L_mat %>% SizeFactorCal()

cellCounts <- read.table("../data/cellCounts.txt", header = T)

cellCounts$L_sf <- tx_LRNA_sf[match(cellCounts$Samples, names(tx_LRNA_sf))] %>% log()
cellCounts$F_sf <- tx_FRNA_sf[match(cellCounts$Samples, names(tx_FRNA_sf))] %>% log()

cellCounts$Color <- factor(gsub("(.*)\\_rep.*", "\\1", cellCounts$Samples),
                           c("SL", "2i_2d", "SL2i_2d", "2i_7d", "mTORi_1d", "mTORi_2d"))

ggplot(cellCounts, aes(x = CellNum / 145, y = L_sf - F_sf, color = Color)) + 
  geom_rect(aes(xmin=0, xmax=0.14, ymin=-1.2, ymax=1), 
            fill="grey80", alpha = 0.3, linetype = 0) +
  geom_line(cex = 1.4, lty = 2) +
  geom_point(cex = 5) + 
  geom_text(aes(x = CellNum / 145, y = L_sf - F_sf, 
                label = gsub(".*_rep(.)", "\\1", Samples)),
            nudge_x = 0.015, color = "black", size = 5) +
  ylim(c(-1.2, 1)) +
  xlab(expression(paste("Cell density (million / ", cm^2, ")"))) +
  ylab("Log cell level turnover") +
  scale_colour_manual(values = colors_20[c(13, 2, 10, 20, 7)]) +
  labs(color='Samples') +
  theme_setting
ggsave(filename = "Fig1.Turn_over_cell_number.pdf", 
       path = "../fig1/figs",
       device = "pdf", width = 6, height = 4 )

tx_FRNA_sf <- tx_F_mat[, !grepl("SL_rep3", colnames(tx_F_mat))] %>% 
  `colnames<-`(gsub("(.*)\\_rep.*","\\1", colnames(.))) %>% 
  meanSampleCounts() %>%
  SizeFactorCal()

tx_LRNA_sf <- tx_L_mat[, !grepl("SL_rep3", colnames(tx_L_mat))] %>%
  `colnames<-`(gsub("(.*)\\_rep.*","\\1", colnames(.))) %>% 
  meanSampleCounts() %>%
  SizeFactorCal()

sf_dat <- data.frame(Samples = names(tx_FRNA_sf),
                     L_sf = log(tx_LRNA_sf),
                     F_sf = log(tx_FRNA_sf))

ggplot(sf_dat, aes(x = F_sf, y = L_sf, color = Samples, label = Samples)) + 
  geom_point(size = 3) +
  geom_text(cex = 5, fontface = "bold", hjust = 0, nudge_x = 0.05) + 
  xlab("Log Total RNA Sample Size") +
  ylab("Log Labeled RNA Sample Size") +
  xlim(c(-1, 1)) + 
  geom_hline(yintercept=0, cex = 0.2) + 
  geom_vline(xintercept=0, cex = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype="dotted") +
  geom_abline(intercept = 0, slope = -1, linetype="dotted") +
  scale_colour_manual(values = colors_20[c(2, 1, 20, 7, 13)]) +
  theme_setting +
  theme(legend.position = "none")
ggsave(filename = "Fig1.size_factors_distribution.pdf", path = "figs",
       device = "pdf", width = 4, height = 4 )


# ------------------------------------------------------------------ #
pca_L <- prcomp(t(log1p(tx_L_mat[names(gene.gr), ])))
pca_F <- prcomp(t(log1p(tx_F_mat[names(gene.gr), ])))
pca_F$x[, "PC2"] <- pca_F$x[, "PC2"] * (-1)

plot_pca <- function(pca, sample_names) {
  pca$sdev <- pca$sdev^2
  pc_var <- pca$sdev[1:2] / sum(pca$sdev) * 100
  pc <- pca$x[, 1:2] %>% as.data.frame()
  pc$Sample <- gsub("_.$", "", sample_names)
  pc$Label <- paste0("rep", gsub(".*_(.$)", "\\1", sample_names))
  pc$Group <- gsub("_.d", "", pc$Sample)
  pc$Group[pc$Group == "SL2i"] <- "2i"
  
  hull <- pc %>%
    group_by(Group) %>%
    slice(chull(PC1, PC2))
  
  ggplot(pc, aes(x = PC1, y = PC2, color = Sample, group = Group)) +
    geom_point(size = 4) +
    geom_text(aes(label = Label), nudge_y = -diff(range(pc$PC2)) / 25) +
    # geom_polygon(data = hull, alpha = 0.5, ) +
    # ggforce::geom_mark_ellipse() +
    # stat_ellipse(aes(x = PC1, y = PC2, group = Group), type = "norm") +
    scale_color_manual(values = colors_20[c(2, 16, 20, 7, 13, 10)]) +
    xlab(paste0("PC1 ", round(pc_var[1]), "% variance")) +
    ylab(paste0("PC2 ", round(pc_var[2]), "% variance")) +
    theme_setting + 
    theme(legend.position = "none")
}

g1 <- plot_pca(pca = pca_L, sample_names = gsub("rep", "", colnames(tx_F_mat))) + ggtitle("Labeled mRNA (10674)")
# g2 <- plot_pca(pca = pca_L_tu, sample_names = gsub("rep", "", colnames(tx_F_mat))) + ggtitle("Labeled ncRNA (41870)")
g3 <- plot_pca(pca = pca_F, sample_names = gsub("rep", "", colnames(tx_F_mat))) + ggtitle("Total mRNA (10674)")

ggsave(plot = grid.arrange(g1, g3, nrow = 1),
       filename = "Fig1.PCA_labeled_total_mRNA.png", 
       path = "../fig1/figs/",
       device = "png", width = 8, height = 4 )

g3 <- plot_pca(pca = pca_F, sample_names = gsub("rep", "", colnames(tx_F_mat))) + ggtitle("Total mRNA (10674)")
g4 <- plot_pca(pca = pca_F_tu, sample_names = gsub("rep", "", colnames(tx_F_mat))) + ggtitle("Total ncRNA (41870)")

ggsave(plot = grid.arrange(g3, g4, nrow = 1),
       filename = "Fig1.PCA_total_mRNA_ncRNA.png", 
       path = "../fig1/figs/",
       device = "png", width = 8, height = 4 )
