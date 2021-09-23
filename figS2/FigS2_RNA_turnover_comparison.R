# Compare samples RNA turnover, by transcript types, gene ontology
# Rui Shao Jun 2020

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(gridExtra)
library(grid)

# gene name convertion table -------------------------------------------------------------------------------
res_biotypes <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                keys = as.character(unique(sample_Tx_counts_Rates$gene_id)),
                                keytype = "GENEID",
                                columns = c("GENENAME", "GENEBIOTYPE", "TXBIOTYPE"))

# Scatter plots of copy number ~ half-life --------------------------------------------------------------------
# GENCODE transcripts
sample_Tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined2.RData")
dat <- sample_Tx_counts_Rates %>%
  dplyr::filter(Sample == "SL") %>%
  dplyr::filter(gene_id %in% res_biotypes$GENEID[res_biotypes$GENEBIOTYPE == "protein_coding"]) %>%
  dplyr::filter(FRNA > 0 & LRNA > 0) %>% 
  dplyr::filter(!is.infinite(Copy) & !is.na(Copy)& Copy > 0) %>%
  dplyr::filter(!is.infinite(Half_life) & !is.na(Half_life) & Half_life < 2000 & Half_life > 1)
g0.1 <- ggplot(dat, aes(x = log10(Copy), y = log10(Half_life))) +
  # geom_point(aes(color  = get_dens(log10(Copy), log10(Half_life)))) +
  geom_hex(binwidth=c(0.1, 0.05)) +
  geom_text(x = 1.8, y = 3.5, label = paste("r =", round(cor(log10(dat$Half_life), log10(dat$Copy)),3)), cex = 4.5, hjust = 0) +
  geom_text(x = 1.8, y = 3.2, label = paste("n =", nrow(dat) ), cex = 4.5, hjust = 0) +
  geom_text(x = -1.6, y = 3.5, label = "mRNA", cex = 4.5, hjust = 0) +
  geom_vline(xintercept = median(log10(dat$Copy), na.rm = T), color = "grey", lty = 2) + # 7.80945
  geom_hline(yintercept = median(log10(dat$Half_life), na.rm = T), color = "grey", lty = 2) + # 56.45986
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="Copy (per cell)", breaks=c(-1, 0, 1, 2), labels=c(0.1, 1, 10, 100), limits = c(-1.5, 3.2)) +
  scale_y_continuous(name="Turnover (min)", breaks=0:3, labels=10^(0:3), limits = c(0, 3.5)) +
  labs(fill = "Genes", title = "Serum LIF") +
  theme(legend.position = "none")

# TU annotations
sample_ncRNA_counts_Rates <- readRDS("../figS2/data/sample_nc_counts_Rates_combined2.RData")
sample_ncRNA_counts_Rates$location <- gsub("_.*", "\\1", sample_ncRNA_counts_Rates$gene_id)
sample_ncRNA_counts_Rates$location[sample_ncRNA_counts_Rates$location == "antisense"] <- "asRNA"

plot_scatter_copy_half_life <- function(dat, title, color_option, .location) {
  print(paste("Median half-life:", median(dat$Half_life, na.rm = T)))
  print(paste("Median copy:", median(dat$Copy, na.rm = T)))
  
  ggplot(dat, aes(x = log10(Copy), y = log10(Half_life))) +
    geom_hex(binwidth=c(0.1, 0.05)) +
    geom_text(x = 1.8, y = 3.5, label = paste("r =", round(cor(log10(dat$Half_life), log10(dat$Copy)),3) ), cex = 4.5, hjust = 0) +
    geom_text(x = -1.6, y = 3.5, label = .location, cex = 4.5, hjust = 0) +
    geom_text(x = 1.8, y = 3.2, label = paste("n =", nrow(dat)), cex = 4.5, hjust = 0) +
    geom_vline(xintercept = median(log10(dat$Copy), na.rm = T), color = "grey", lty = 2) +
    geom_hline(yintercept = median(log10(dat$Half_life), na.rm = T), color = "grey", lty = 2) +
    theme_setting +
    scale_fill_viridis(option = color_option, direction = -1) +
    scale_x_continuous(name="Copy (per cell)", breaks=c(-1:2), labels=10^c(-1:2), limits = c(-1.5, 3)) +
    scale_y_continuous(name="Half-life (min)", breaks=0:3, labels=10^(0:3), limits = c(0, 3.5)) +
    labs(fill = "TUs", title = title) +
    theme(legend.position = "none")
}

make_ggplot_object <- function(table, title, .Sample, .location, color_option = "A") {
  table %>%
    dplyr::filter(Sample == .Sample) %>%
    dplyr::filter(grepl(.location, location)) %>%
    dplyr::filter(!is.infinite(Copy) & !is.na(Copy)& Copy > 0) %>%
    dplyr::filter(!is.infinite(Half_life) & !is.na(Half_life) & Half_life < 2000 & Half_life > 1) %>%
    plot_scatter_copy_half_life(title = title, color_option, .location)
}

g0.2 <- make_ggplot_object(sample_ncRNA_counts_Rates,
                           title = "",
                           .Sample = "SL",
                           .location = "intergenic") # "Median half-life: 11.4218252647328" "Median copy: 2.08828830408298"
g0.3 <- make_ggplot_object(sample_ncRNA_counts_Rates,
                           title = "",
                           .Sample = "SL",
                           .location = "uaRNA") # "Median half-life: 7.69946326879276" "Median copy: 1.2210601025218"
g0.4 <- make_ggplot_object(sample_ncRNA_counts_Rates,
                           title = "",
                           .Sample = "SL",
                           .location = "asRNA") # "Median half-life: 6.60008578309882" "Median copy: 2.61632784323352"

ggsave(plot = grid.arrange(g0.1, g0.2, g0.3, g0.4, nrow = 2),
       filename = "FigS2_Scatter_halflife_copy_SL.png",
       path = "../figS2/figs", device = "png", width = 8, height = 8)

# mRNA condition comparisons ----------------------------------------------------------------------
copy_table <- data.frame("SL" = sample_Tx_counts_Rates$Copy[sample_Tx_counts_Rates$Sample == "SL"],
                         "2i" = sample_Tx_counts_Rates$Copy[sample_Tx_counts_Rates$Sample == "2i_2d"],
                         "mTORi" = sample_Tx_counts_Rates$Copy[sample_Tx_counts_Rates$Sample == "mTORi_1d"])
half_life_table <- data.frame("SL" = sample_Tx_counts_Rates$Half_life[sample_Tx_counts_Rates$Sample == "SL"],
                              "2i" = sample_Tx_counts_Rates$Half_life[sample_Tx_counts_Rates$Sample == "2i_2d"],
                              "mTORi" = sample_Tx_counts_Rates$Half_life[sample_Tx_counts_Rates$Sample == "mTORi_1d"])
# keep mRNA and remove NAs
copy_table <- copy_table[unique(sample_Tx_counts_Rates$gene_id) %in% res_biotypes$GENEID[res_biotypes$GENEBIOTYPE == "protein_coding"], ]
copy_table <- copy_table[!apply(copy_table, 1, function(x) any(is.na(x) | is.infinite(x))), ]
half_life_table <- half_life_table[unique(sample_Tx_counts_Rates$gene_id) %in% res_biotypes$GENEID[res_biotypes$GENEBIOTYPE == "protein_coding"], ]
half_life_table <- half_life_table[!apply(half_life_table, 1, function(x) any(is.na(x) | is.infinite(x))), ]


g1.1 <- ggplot(copy_table, aes(y = log10(SL), x = log10(X2i))) +
  geom_hex(bins = 40) +
  geom_abline(intercept = 0, color = "grey", lty = 2) +
  geom_text(x = -2, y = 4, label = paste("ΔMedian =", round(median(copy_table$X2i - copy_table$SL), 2) ), cex = 4.5, hjust = 0) +
  geom_text(x = -2, y = 3.6, label = paste("n =", nrow(copy_table) ), cex = 4.5, hjust = 0) +
  theme_setting +
  scale_fill_viridis(option = "E", direction = -1, end = 0.9) +
  scale_y_continuous(name="Serum LIF", breaks=c(-1:3), labels=10^(-1:3), limits = c(-2, 4)) +
  scale_x_continuous(name="2i 2d", breaks=(-1:3), labels=10^(-1:3), limits = c(-2, 4)) +
  labs(title = "Copy") +
  theme(legend.position = "none")

g1.2 <- ggplot(copy_table, aes(y = log10(SL), x = log10(mTORi))) +
  geom_hex(bins = 40) +
  geom_abline(intercept = 0, color = "grey", lty = 2) +
  geom_text(x = -2, y = 4, label = paste("ΔMedian =", round(median(copy_table$mTORi - copy_table$SL), 2) ), cex = 4.5, hjust = 0) +
  geom_text(x = -2, y = 3.6, label = paste("n =", nrow(copy_table) ), cex = 4.5, hjust = 0) +
  theme_setting +
  scale_fill_viridis(option = "D", direction = -1, end = 0.85) +
  scale_y_continuous(name="Serum LIF", breaks=c(-1:3), labels=10^(-1:3), limits = c(-2, 4)) +
  scale_x_continuous(name="mTORi 1d", breaks=(-1:3), labels=10^(-1:3), limits = c(-2, 4)) +
  labs(title = "Copy") +
  theme(legend.position = "none")

g1.3 <- ggplot(half_life_table, aes(y = log10(SL), x = log10(X2i))) +
  geom_hex(bins = 40) +
  geom_abline(intercept = 0, color = "grey", lty = 2) +
  geom_text(x = 0, y = 3.5, label = paste("ΔMedian =", round(median(half_life_table$X2i - half_life_table$SL), 2) ), cex = 4.5, hjust = 0) +
  geom_text(x = 0, y = 3.25, label = paste("n =", nrow(half_life_table) ), cex = 4.5, hjust = 0) +
  theme_setting +
  scale_fill_viridis(option = "E", direction = -1, end = 0.9) +
  scale_y_continuous(name="Serum LIF", breaks=c(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  scale_x_continuous(name="2i 2d", breaks=(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  labs(title = "Half-life (min)") +
  theme(legend.position = "none")

g1.4 <- ggplot(half_life_table, aes(y = log10(SL), x = log10(mTORi))) +
  geom_hex(bins = 40) +
  geom_abline(intercept = 0, color = "grey", lty = 2) +
  geom_text(x = 0, y = 3.5, label = paste("ΔMedian =", round(median(half_life_table$mTORi - half_life_table$SL), 2) ), cex = 4.5, hjust = 0) +
  geom_text(x = 0, y = 3.25, label = paste("n =", nrow(half_life_table) ), cex = 4.5, hjust = 0) +
  theme_setting +
  scale_fill_viridis(option = "D", direction = -1, end = 0.85) +
  scale_y_continuous(name="Serum LIF", breaks=c(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  scale_x_continuous(name="mTORi 1d", breaks=(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  labs(title = "Half-life (min)") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1.1, g1.3, g1.2, g1.4, nrow = 1),
       filename = "FigS2_halflife_copy_comparisons.png", path = "figs",
       device = "png", width = 16, height = 4)

# compare different TU types between conditions ------------------------------------------------------------------------
sample_Tx_counts_Rates$location <- "mRNA"
dat <- rbind(sample_Tx_counts_Rates[sample_Tx_counts_Rates$FRNA > 0, ],
             sample_ncRNA_counts_Rates)
dat$Sample <- factor(dat$Sample, c("SL", "2i_2d", "2i_7d", "mTORi_1d", "mTORi_2d"))
dat$location <- factor(dat$location, c("mRNA", "intergenic", "uaRNA", "asRNA", "conRNA", "daRNA", "usRNA", "dsRNA", "gene"))
dat <- dat %>% dplyr::filter(!is.na(Copy) & !is.infinite(Copy) & Copy > 0.01 &
                               !is.na(Half_life) & !is.infinite(Half_life))


g1.5 <- ggplot(dat[dat$location %in% c("mRNA", "intergenic", "uaRNA", "asRNA", "conRNA") &
             dat$Sample %in% c("SL", "2i_2d", "mTORi_1d"), ], 
       aes(x = location, y = log10(Half_life), fill = Sample)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = colors_20[c(13, 2, 20)]) +
  xlab("") +
  scale_y_continuous(name="Turnover half-life (min)", breaks=0:3, labels=10^(0:3), limits = c(0, 3.2)) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")

g1.6 <- ggplot(dat[dat$location %in% c("mRNA", "intergenic", "uaRNA", "asRNA", "conRNA") &
               dat$Sample %in% c("SL", "2i_2d", "mTORi_1d"), ], 
       aes(x = location, y = log10(Copy), fill = Sample)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = colors_20[c(13, 2, 20)]) +
  xlab("") + 
  scale_y_continuous(name="Copy (per cell)", breaks=(-1):2, labels=10^((-1):2), limits = c(-1.2, 2.7)) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(plot = grid.arrange(g1.5, g1.6, nrow = 1, widths = c(3, 4.4)),
       filename = "FigS2_sample_copy_half_life_boxplot.png", path = "../figS2/figs",
       device = "png", width = 7, height = 4)

# Plot log2FC comparison --------------------------------------------------------------------
TU.DE.mm10 <- readRDS("../fig1/data/TU.DE.mm10.RData")
TU.DE.mm10 <- TU.DE.mm10[apply(TU.DE.mm10[, -1], 1, function(x) all(!is.na(x) & is.finite(x))), ]

g2.1 <- ggplot(TU.DE.mm10, aes(x = L_log2FC_mTORi, y = L_log2FC_2i, 
                      color = get_dens(L_log2FC_mTORi, L_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TU.DE.mm10) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Nascent RNA") +
  theme(legend.position = "none")


g2.2 <- ggplot(TU.DE.mm10, aes(x = F_log2FC_mTORi, y = F_log2FC_2i, 
                              color = get_dens(F_log2FC_mTORi, F_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TU.DE.mm10) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Total RNA") +
  theme(legend.position = "none")

g2.3 <- ggplot(TU.DE.mm10, aes(x = L_log2FC_mTORi - F_log2FC_mTORi, 
                      y = L_log2FC_2i - F_log2FC_2i,
                      color = get_dens(L_log2FC_mTORi - F_log2FC_mTORi, L_log2FC_2i - F_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TU.DE.mm10) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC contrast mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC contrast 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Nascent RNA - Total RNA") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, nrow = 1),
       filename = "FigS2_LRNA_FRNA_log2FC_comparisons.png", path = "figs",
       device = "png", width = 12, height = 4)

# Plot log2FC vs half-life changes --------------------------------------------------------------------

# copy change ~ half-life change ------------------------------------------------------------------------
# pluripotent gene list refers to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5399622/
# and https://genome.cshlp.org/content/early/2018/08/28/gr.233437.117.full.pdf

pluripotent_markers <- list("General pluripotency" = c("Myc", "Mycn", "Nanog", "Pou5f1", "Sox2", "Sall4", "Tdgf1", "Utf1", "Zic2", "Zic3"),
                            "Naïve pluripotency" = c("Dnmt3l", "Esrrb", "Klf2", "Klf4", "Klf5", "Nr5a2", "Nr0b1", 
                                                      "Tead4", "Tfcp2l1", "Prdm14", "Zfp42", "Gbx2", "Tbx3"),
                            "Post-implantation" = c("Cdh2", "Dnmt3a", "Dnmt3b", "Etv4", "Etv5", "Fgf5", "Fgf10", "Fgf15", "Fgfr1", "Fzd7",  
                                                    "Id1", "Lefty1", "Lefty2", "Otx2", "Sall2", "Sox3", "Sox4", "Sox11", "Pou3f1", "Nodal"),
                            "Lineage markers" = c("Kit", "Dazl", "Gdf9", "Sox1", "Pax6", "Neurog1", "Neurog2", "Ascl1", "Zic1", "Eomes",
                                                  "Gata4", "Gata6", "Cdx2", "Meox1", "Foxf1", "Foxa2", "Sox7", "Sox17") )

res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(unique(sample_Tx_counts_Rates$gene_id)),
                       keytype = "GENEID",
                       columns = "SYMBOL")
rownames(res) <- res$GENEID

sample_Tx_counts_Rates$gene_symbol <- res[match(sample_Tx_counts_Rates$gene_id, res$GENEID), 2]

dat_markers <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$Sample %in% c("SL", "2i_2d", "mTORi_1d") & 
                                        sample_Tx_counts_Rates$gene_symbol %in% unlist(pluripotent_markers), ]
dat_markers$Type <- rep(names(pluripotent_markers), lengths(pluripotent_markers))[match(dat_markers$gene_symbol, unlist(pluripotent_markers))]

dat_markers <- dat_markers[!is.infinite(dat_markers$Copy), ]
dat_markers <- dat_markers[dat_markers$gene_symbol %in%
                             names(table(dat_markers$gene_symbol)[table(dat_markers$gene_symbol) == 3]), ]
dat_markers$gene_symbol <- factor(dat_markers$gene_symbol,
                                  levels = unique(dat_markers$gene_symbol)[(order(dat_markers[dat_markers$Sample == "SL", "Copy"]))] )
dat_markers$Sample <- factor(dat_markers$Sample, c("SL", "2i_2d", "mTORi_1d"))
dat_markers$Type <- factor(dat_markers$Type, names(pluripotent_markers))

# by dot charts
g3.1 <- ggplot(dat_markers[dat_markers$Type %in% c("General pluripotency", "Naïve pluripotency"), ], 
             aes(x = gene_symbol, y = log10(Copy), 
                        color = Sample, group = gene_symbol)) +
  geom_line(color = "grey50") +
  geom_point(cex = 2, shape=21, stroke = 1.5) +
  ylab("Copy / cell") +
  xlab("") +
  scale_y_continuous(breaks = (0:3), labels = 10^(0:3), limits = c(-2, 3)) +
  scale_color_manual(values = colors_20[c(13,2,20)]) +
  facet_grid(.~ Type, switch = "y", scales = "free", space = "free") + 
  theme_setting + 
  theme_minimal() +
  theme(strip.placement = "outside",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=11, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_line(),
        strip.background = element_rect(), 
        strip.text.x = element_text(size=11),
        legend.position = "none")

g3.2 <- ggplot(dat_markers[!dat_markers$Type %in% c("General pluripotency", "Naïve pluripotency"), ], 
             aes(x = gene_symbol, y = log10(Copy), 
                 color = Sample, group = gene_symbol)) +
  geom_line(color = "grey50") +
  geom_point(cex = 2, shape=21, stroke = 1.5) +
  ylab("Copy / cell") +
  xlab("") +
  scale_y_continuous(breaks = (-1:2), labels = 10^(-1:2), limits = c(-2, 2)) +
  scale_color_manual(values = colors_20[c(13,2,20)]) +
  facet_grid(.~ Type, switch = "y", scales = "free", space = "free") + 
  theme_setting + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=11, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_line(),
        strip.placement = "outside",
        strip.background = element_rect(), 
        strip.text.x = element_text(size=11),
        legend.position = "none")

ggsave(plot = grid.arrange(g3.1, g3.2, nrow = 2),
       filename = "FigS2_pluripotent_markers_copy_.png",
       path = "../figS2/figs", device = "png", width = 8, height = 5)

# by heatmap
dat_markers$Sample <- factor(dat_markers$Sample, rev(c("SL", "2i_2d", "mTORi_1d")))
g3.1 <- ggplot(dat_markers[dat_markers$Type %in% c("General pluripotency", "Naïve pluripotency"), ], 
       aes(x = gene_symbol, y = Sample, fill = log10(Copy + 1))) +
  geom_tile() +
  facet_grid(.~Type, scales = "free", space = "free") +
  scale_fill_gradientn(colours = rev(brewer.pal(10, "RdBu"))) +
  xlab("") + ylab("") +
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key.size = unit(15, 'pt'))
  
g3.2 <- ggplot(dat_markers[!dat_markers$Type %in% c("General pluripotency", "Naïve pluripotency"), ], 
               aes(x = gene_symbol, y = Sample, fill = log10(Copy + 1))) +
  geom_tile() +
  facet_grid(.~Type, scales = "free", space = "free") +
  scale_fill_gradientn(colours = rev(brewer.pal(10, "RdBu"))) +
  xlab("") + ylab("") +
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.key.size = unit(15, 'pt'))

ggsave(plot = grid.arrange(g3.1, g3.2, nrow = 2),
       filename = "FigS2_heatmap_pluripotent_markers_copy.png",
       path = "../figS2/figs", device = "png", width = 10, height = 6)


# contribution of copy change, half-life
h.idx <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$Sample == "SL", "FRNA"] > 50
l.idx <- sample_Tx_counts_Rates[sample_Tx_counts_Rates$Sample == "SL", "FRNA"] < 1

ratio_changes <- function(X, Y)
  ifelse((X==0 | Y==0 | (X*Y < 0.001)), NA, 1) * (log2(X + 1) - log2(Y + 1)) / ((log2(X + 1) + log2(Y + 1)) / 2)

log2FC_cmp <- function(X1, Y1, X2, Y2)
  ifelse(X1 * Y1 * X2 * Y2 < 0.001, NA, 1) * (log2(X1 + 1) - log2(Y1 + 1)) - ((log2(X2 + 1) - log2(Y2 + 1)))


rho <- readRDS("../fig2/data/est.rho.RData")

dat_rho <- with(sample_Tx_counts_Rates[sample_Tx_counts_Rates$Sample %in% c("SL", "2i_2d", "mTORi_1d"), ],
                cbind(matrix(Copy, nrow = sum(sample_Tx_counts_Rates$Sample == "SL")),
                      matrix(Labeled_rate, nrow = sum(sample_Tx_counts_Rates$Sample == "SL")),
                      matrix(Half_life, nrow = sum(sample_Tx_counts_Rates$Sample == "SL")),
                      matrix(FRNA, nrow = sum(sample_Tx_counts_Rates$Sample == "SL")),
                      matrix(LRNA, nrow = sum(sample_Tx_counts_Rates$Sample == "SL")) ) )
colnames(dat_rho) <- c(sapply(c("Copy", "Labeled_rate", "Half_life", "FRNA", "LRNA"), 
                            function(x) paste(x, c("2i_2d", "mTORi_1d", "SL"), sep = "_")))
rownames(dat_rho) <- unique(sample_Tx_counts_Rates$gene_id)
dat_rho <- dat_rho[names(rho), ] %>% as.data.frame()
dat_rho$rho <- rho

g3.1 <- ggplot(dat_rho, aes(x = rho, y = log10(Copy_SL), 
                    color = get_dens(rho, log10(Copy_SL), 300))) +
  geom_point(size = 1) + 
  geom_smooth(method = "loess", lwd = 0.5,  color = "red2", se = F) +
  scale_color_gradient(low = add.alpha("grey80", 0.5), high = "grey20") +
  xlim(c(0, 1)) +
  xlab("Saturation ρ") +
  ylab("log10 Copy (cell)") +
  theme_setting +
  theme(legend.position = "none")

g3.2 <- ggplot(dat_rho, aes(x = rho, y = log10(Half_life_SL), 
                    color = get_dens(rho, log10(Half_life_SL), 300))) +
  geom_point(size = 1) + 
  geom_smooth(method = "loess", lwd = 0.5,  color = "red2", se = F) +
  scale_color_gradient(low = add.alpha("grey80", 0.5), high = "grey20") +
  xlim(c(0, 1)) +
  xlab("Saturation ρ") +
  ylab("log10 Turnover half-life (min)") +
  theme_setting +
  theme(legend.position = "none")

g3.3 <- ggplot(dat_rho, aes(x = rho, y = log10(Labeled_rate_SL / 5 * Copy_SL), 
                    color = get_dens(rho, log10(Labeled_rate_SL / 5 * Copy_SL), 300))) +
  geom_point(size = 1) + 
  geom_smooth(method = "loess", lwd = 0.5,  color = "red2", se = F) +
  scale_color_gradient(low = add.alpha("grey80", 0.5), high = "grey20") +
  xlim(c(0, 1)) +
  xlab("Saturation ρ") +
  ylab("log10 Synthesis rate (copy / min)") +
  theme_setting +
  theme(legend.position = "none")

margin = theme(plot.margin = unit(c(0,1,0,1), "cm"))

ggsave(plot = grid.arrange(grobs = lapply(list(g3.2, g3.1, g3.3), "+", margin),nrow = 1),
       filename = "FigS2_saturation_copy_half_life_mu.png",
       path = "../figS2/figs", device = "png", width = 12, height = 3.5)

modes <- function(d) {
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

# plot rho by pluripotent gene group
dat_rho <- data.frame("rho" = rho, 
                      "Type" = dat_markers$Type[match(names(rho), dat_markers$gene_id)],
                      "gene_name" = dat_markers$gene_symbol[match(names(rho), dat_markers$gene_id)]) %>% 
            dplyr::filter(complete.cases(.))

x_axis_labs <- c("General\npluripotency", "Naïve\npluripotency", "Post\nimplantation", "Lineage\nmarkers")

g4.1 <- ggplot(dat_markers, aes(x = Type, y = log10(Copy), color = Type)) +
  geom_boxplot(outlier.alpha = 0, lwd = 1) +
  geom_jitter(width = 0.1, size = 1) + 
  xlab("") + ylab( expression(paste("Copy (cell)"^"-1")) ) +
  scale_color_viridis_d(end = 0.7) +
  scale_x_discrete(labels = x_axis_labs) +
  scale_y_continuous(breaks = (-1:2), labels = 10^(-1:2), limits = c(-1, 2.8)) +
  theme_setting +
  theme(legend.position = "none",
        axis.text.x=element_blank())

g4.2 <- ggplot(dat_markers, aes(x = Type, y = log10(Half_life), color = Type)) +
  geom_boxplot(outlier.alpha = 0, lwd = 1) +
  geom_jitter(width = 0.1, size = 1) + 
  xlab("") + ylab("Turnover half-life (min)") +
  scale_color_viridis_d(end = 0.7) +
  scale_x_discrete(labels = x_axis_labs) +
  scale_y_continuous(breaks = (0:3), labels = 10^(0:3), limits = c(0.7, 3)) +
  theme_setting +
  theme(legend.position = "none",
        axis.text.x=element_blank())

g4.3 <- ggplot(dat_markers, aes(x = Type, y = log10(Labeled_rate * Copy), color = Type)) +
  geom_boxplot(outlier.alpha = 0, lwd = 1) +
  geom_jitter(width = 0.1, size = 1) + 
  xlab("") + ylab( expression(paste("Synthesis (cell·min)"^"-1")) ) +
  scale_color_viridis_d(end = 0.7) +
  scale_x_discrete(labels = x_axis_labs) +
  scale_y_continuous(breaks = (-3:2), labels = 10^(-3:2), limits = c(-2.4, 1.7)) +
  theme_setting +
  theme(legend.position = "none",
        axis.text.x=element_blank())

g4.4 <- ggplot(dat_rho, aes(x = Type, y = log10(rho), color = Type)) +
  geom_boxplot(outlier.alpha = 0, lwd = 1) +
  geom_jitter(width = 0.1, size = 1) + 
  xlab("") + ylab("Saturation ρ") +
  scale_color_viridis_d(end = 0.7) +
  scale_x_discrete(labels = x_axis_labs) +
  scale_y_continuous(breaks = c(-1, -0.5, 0), 
                     labels = round(10^c(-1, -0.5, 0), 1),
                     limits = c(-1.2, 0.4)) +
  theme_setting +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(plot = cowplot::plot_grid(g4.1, g4.2, g4.3, g4.4, ncol = 1, align = "v"),
       filename = "FigS2_pluripotent_group_copy_half_life_mu.png",
       path = "../figS2/figs", device = "png", width = 3.5, height = 10)

# ---------------------------------------------------------------------------------------------------
# compare m6A profile

filenames <- sort(list.files('../data/kallisto_output_m6A/', ignore.case = T, full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "est_counts")
colnames(count_table) <- sampleNewName
count_table <- count_table[, grep("Abcam", colnames(count_table))]
count_table <- rbind(keepOneTx(count_table, rowname_gene_id = T),
                     count_table[!grepl("^chrS|^ENS", rownames(count_table)), ])

library(DESeq2)
dds_m6A <- DESeqDataSetFromMatrix(round(as.matrix(count_table)),
                                   colData = data.frame(condition = gsub(".*Abcam_(.*)_Rep.", "\\1", colnames(count_table)) ),
                                   design = ~ condition)
dds_m6A <- DESeq(dds_m6A)
res_m6A <- results(dds_m6A, contrast = c("condition", "ChrMeRIP", "Input"))

dat_TU_rate <- rbind(sample_Tx_counts_Rates[sample_Tx_counts_Rates$Sample == "SL", ],
                sample_ncRNA_counts_Rates[sample_ncRNA_counts_Rates$Sample == "SL", ])
dat_TU_rate <- dat_TU_rate[match(rownames(res_m6A), dat_TU_rate$gene_id), ]
dat_TU_rate$m6A <- res_m6A$log2FoldChange

dat_TU_rate <- dat_TU_rate %>% 
  filter( Copy > 1 & location %in% c("mRNA", "intergenic") ) %>%
  filter(!is.infinite(Half_life) & !is.na(Half_life))

dat_TU_rate$Half_life_grid <- cut(dat_TU_rate$Half_life, 
                                  c(0, 10, 20, 30, 60, Inf))
# t.test mRNA vs intergenic m6A: 
#   (30,60], p-value = 0.02168
#   (60,Inf], p-value = 9.148e-05

ggplot(dat_TU_rate, aes(x = Half_life_grid, y = m6A, fill = location)) +
  geom_hline(yintercept = 0, color = "grey50", size = 1) +
  geom_violin(position = position_dodge(width = 0.6))+
  geom_boxplot(outlier.size = 0, width=.4, position = position_dodge(width = 0.6)) +
  annotate("text",x = 1, y = 5, label = "ns") +
  annotate("text",x = 2, y = 5, label = "ns") +
  annotate("text",x = 3, y = 5, label = "ns") +
  annotate("text",x = 4, y = 5, label = "*") +
  annotate("text",x = 5, y = 5, label = "****") +
  scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.9, alpha = 0.8) +
  ylim(c(-5, 5)) +
  xlab("\nTurnover half-life (min)") + ylab("m6A enrichment\n") +
  theme_setting

ggsave(filename = "FigS2_mRNA_intergenic_turnover_m6A.png",
       device = "png", path = "../figS2/figs", width = 6, height = 4)

# compare across enhancer types
intergenic.mm10.gr <- rownames(count_table)[grepl("intergenic", rownames(count_table))] %>%
  strsplit("_") %>%
  Reduce("rbind", x = .) %>%
  as.data.frame() %>%
  `colnames<-`(c("location", "seqname", "start", "end", "strand")) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

intergenic.mm10.gr$gene_id <- rownames(count_table)[grepl("intergenic", rownames(count_table))]
intergenic.mm10.gr$enhancer_type <- "none"
mtch <- findOverlaps(enhancer_mm10, intergenic.mm10.gr)
intergenic.mm10.gr$enhancer_type[subjectHits(mtch)] <- enhancer_mm10$type[queryHits(mtch)]

dat_intergenic_rate <- sample_ncRNA_counts_Rates %>% dplyr::filter(Sample == "SL" & location == "intergenic")
dat_intergenic_rate <- dat_intergenic_rate[match(intergenic.mm10.gr$gene_id, dat_intergenic_rate$gene_id), ]

dat_intergenic_rate$m6A <- res_m6A$log2FoldChange[grepl("intergenic", rownames(count_table))]
dat_intergenic_rate$Enhancer <- factor(intergenic.mm10.gr$enhancer_type, levels = c("active", "pseudo", "clean", "none"))
dat_intergenic_rate <- dat_intergenic_rate %>% filter(complete.cases(.) & !is.infinite(Half_life))

stat.test <- dat_intergenic_rate %>%
  mutate(Half_life = log10(Half_life)) %>%
  t_test(Half_life ~ Enhancer) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(x = "Enhancer")
stat.test <- stat.test[1:3, ]

library(ggpubr)
library(rstatix)
g5.1 <- ggboxplot(dat_intergenic_rate, x = "Enhancer", y = "Half_life", fill = "Enhancer", outlier.shape = NA) +
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "1e1", "1e2", "1e3"), limits = c(1, 1000)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = c(2.6, 2.8, 3), tip.length = 0.01) +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 1, alpha = 0.8) +
  xlab("") + ylab("Turnover half-life (min)\n") +
  theme_setting + 
  theme(axis.text.x = element_blank())

stat.test <- dat_intergenic_rate %>%
  t_test(m6A ~ Enhancer) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(x = "Enhancer")
stat.test <- stat.test[1:3, ]

g5.2 <- ggboxplot(dat_intergenic_rate, x = "Enhancer", y = "m6A", fill = "Enhancer", outlier.shape = NA) +
  ylim(c(-3, 5)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = c(3.8, 4.4, 5), tip.length = 0.01) +
  scale_fill_viridis_d(option = "B", begin = 0.3, end = 1, alpha = 0.8) +
  xlab("") + ylab("m6A enrichment\n") +
  theme_setting +
  theme(axis.text.x = element_blank())

ggsave(plot = grid.arrange(g5.1, g5.2, ncol=1),
       filename = "FigS2_intergenic_enhancer_type_m6A_turnover_m6A.png",
       device = "png", path = "../figS2/figs", width = 5, height = 6)


# compare to TT-seq 2016 science paper ------------------------------------------------------------------------
k562_anno_hg38.gr <- importRanges("../data/GSE75792_transcript.annotation.gtf")
k562_anno_hg38.gr$synthesis <- gsub('rate "', "", k562_anno_hg38.gr$synthesis) %>% as.numeric()
k562_anno_hg38.gr$decay <- gsub('rate "', "", k562_anno_hg38.gr$decay) %>% as.numeric()

gene.gr <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
gene.gr <- gene.gr[gene.gr$gene_biotype == "protein_coding"]

mtch <- findOverlaps(k562_anno_hg38.gr, gene.gr)
k562_anno_hg38.gr$gene_id <- NA
k562_anno_hg38.gr$gene_id[queryHits(mtch)] <- gene.gr$gene_id[subjectHits(mtch)]

k562_anno_hg38.gr <- k562_anno_hg38.gr[k562_anno_hg38.gr$type == "protein_coding"]

genesV2 = biomaRt::getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                 values = k562_anno_hg38.gr$gene_id,
                 mart = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl"),
                 martL = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                 attributesL = c("ensembl_gene_id"), 
                 uniqueRows = T)
k562_anno_hg38.gr$m_gene_id <- genesV2[, 2][match(k562_anno_hg38.gr$gene_id, genesV2[, 1])]
k562_anno_hg38.gr <- k562_anno_hg38.gr[!is.na(k562_anno_hg38.gr$decay) & !is.na(k562_anno_hg38.gr$m_gene_id)]
k562_anno_hg38.gr$half_life <- log(2) / k562_anno_hg38.gr$decay

sample_Tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined.RData")
Tx_Rates_SL <- sample_Tx_counts_Rates %>%
  dplyr::filter(Sample == "SL") %>%
  dplyr::filter(gene_id %in% res_biotypes$GENEID[res_biotypes$GENEBIOTYPE == "protein_coding"]) %>%
  dplyr::filter(FRNA > 0 & LRNA > 0) %>% 
  dplyr::filter(!is.infinite(Copy) & !is.na(Copy)& Copy > 0) %>%
  dplyr::filter(!is.infinite(Half_life) & !is.na(Half_life) & Half_life < 2000 & Half_life > 1)

data.frame(mESC = log10(Tx_Rates_SL$Half_life[match(k562_anno_hg38.gr$m_gene_id, Tx_Rates_SL$gene_id)]),
     K562 = log10(k562_anno_hg38.gr$half_life)) %>%
  dplyr::filter(complete.cases(.) & !is.infinite(mESC) & !is.infinite(K562)) %>%
  dplyr::mutate("density" = get_dens(mESC, K562, n.grid = 200)) %>%
  ggplot(aes(x = mESC, y = K562, color = density)) +
  geom_point(size = 1) + 
  scale_color_gradient(low = "lightgrey", high = 1) +
  scale_x_continuous(name = "mESC half-life (min)",
                     breaks = c(0:3),
                     labels = 10^(0:3)) +
  scale_y_continuous(name = "K562 half-life (min)",
                     breaks = c(0:4),
                     labels = 10^(0:4), 
                     limits = c(0, 4)) +
  theme_setting +
  theme(legend.position = "none")

ggsave(filename = "FigS2_mESC_K562_half_life_comparison.png", path = "../figS2/figs",
       device = "png", width = 5, height = 4)



