# Compare samples RNA turnover, on transcript types, gene ontology
# Rui Shao Jun 2020

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(gridExtra)
library(grid)

# -------------------------------------------------------------------------------
res_biotypes <- biomaRt::select(EnsDb.Mmusculus.v79,
                                keys = as.character(unique(sample_Tx_counts_Rates$gene_id)),
                                keytype = "GENEID",
                                columns = c("GENENAME", "GENEBIOTYPE", "TXBIOTYPE"))

# Scatter plots of copy number ~ half-life --------------------------------------------------------------------
# GENCODE transcripts
sample_Tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined.RData")
dat <- sample_Tx_counts_Rates %>%
  dplyr::filter(Sample == "SL") %>%
  dplyr::filter(gene_id %in% res_biotypes$GENEID[res_biotypes$GENEBIOTYPE == "protein_coding"]) %>%
  dplyr::filter(FRNA > 0 & LRNA > 0) %>% 
  dplyr::filter(!is.infinite(Copy) & !is.na(Copy)& Copy > 0) %>%
  dplyr::filter(!is.infinite(Half_life) & !is.na(Half_life) & Half_life < 2000 & Half_life > 1)
g0.1 <- ggplot(dat, aes(x = log10(Copy), y = log10(Half_life))) +
  geom_hex(binwidth=c(0.1, 0.05)) +
  geom_text(x = 2, y = 3.5, label = paste("r =", round(cor(log10(dat$Half_life), log10(dat$Copy)),3)), cex = 4.5, hjust = 0) +
  geom_text(x = 2, y = 3.2, label = paste("n =", nrow(dat) ), cex = 4.5, hjust = 0) +
  geom_text(x = -1.6, y = 3.5, label = "mRNA", cex = 4.5, hjust = 0) +
  geom_vline(xintercept = median(log10(dat$Copy), na.rm = T), color = "grey", lty = 2) + # 7.80945
  geom_hline(yintercept = median(log10(dat$Half_life), na.rm = T), color = "grey", lty = 2) + # 56.45986
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="Copy (per cell)", breaks=c(-1, 0, 1, 2), labels=c(0.1, 1, 10, 100), limits = c(-1.5, 3.2)) +
  scale_y_continuous(name="Half life (min)", breaks=0:3, labels=10^(0:3), limits = c(0, 3.5)) +
  labs(fill = "Genes", title = "Serum LIF") +
  theme(legend.position = "none")

# TU annotations
sample_ncRNA_counts_Rates <- readRDS("../figS2/data/sample_nc_counts_Rates_combined.RData")
sample_ncRNA_counts_Rates$location <- gsub("_.*", "\\1", sample_ncRNA_counts_Rates$gene_id)

plot_scatter_copy_half_life <- function(dat, title, color_option, .location) {
  print(paste("Median half-life:", median(dat$Half_life, na.rm = T)))
  print(paste("Median copy:", median(dat$Copy, na.rm = T)))
  
  ggplot(dat, aes(x = log10(Copy), y = log10(Half_life))) +
    geom_hex(binwidth=c(0.1, 0.05)) +
    geom_text(x = 2, y = 2.5, label = paste("r =", round(cor(log10(dat$Half_life), log10(dat$Copy)),3) ), cex = 4.5, hjust = 0) +
    geom_text(x = -1.6, y = 2.5, label = .location, cex = 4.5, hjust = 0) +
    geom_text(x = 2, y = 2.35, label = paste("n =", nrow(dat)), cex = 4.5, hjust = 0) +
    geom_vline(xintercept = median(log10(dat$Copy), na.rm = T), color = "grey", lty = 2) +
    geom_hline(yintercept = median(log10(dat$Half_life), na.rm = T), color = "grey", lty = 2) +
    theme_setting +
    scale_fill_viridis(option = color_option, direction = -1) +
    scale_x_continuous(name="Copy (per cell)", breaks=c(-1:2), labels=10^c(-1:2), limits = c(-1.5, 3)) +
    scale_y_continuous(name="Half life (min)", breaks=0:2, labels=10^(0:2), limits = c(0, 2.5)) +
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

ggsave(plot = grid.arrange(g0.1, g0.2, g0.3, g0.4, nrow = 1),
       filename = "Fig2_Scatter_halflife_copy_SL.png", path = "figs",
       device = "png", width = 16, height = 4)

# condition comparisons
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
  scale_fill_viridis(option = "E", direction = -1) +
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
  scale_fill_viridis(option = "E", direction = -1) +
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
  scale_fill_viridis(option = "E", direction = -1) +
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
  scale_fill_viridis(option = "E", direction = -1) +
  scale_y_continuous(name="Serum LIF", breaks=c(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  scale_x_continuous(name="mTORi 1d", breaks=(0:3), labels=10^(0:3), limits = c(0, 3.5)) +
  labs(title = "Half-life (min)") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1.1, g1.3, g1.2, g1.4, nrow = 1),
       filename = "Fig2_halflife_copy_comparisons.png", path = "figs",
       device = "png", width = 16, height = 4)

# Plot log2FC comparison --------------------------------------------------------------------
TX.DE.tpm <- readRDS("../fig1/data/TX.DE.tpm.external.norm.RData")
TX.DE.tpm <- TX.DE.tpm[apply(TX.DE.tpm[, -1], 1, function(x) all(!is.na(x) & is.finite(x))), ]

g1.6 <- ggplot(TX.DE.tpm, aes(x = L_log2FC_mTORi, y = L_log2FC_2i, 
                      color = get_dens(L_log2FC_mTORi, L_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TX.DE.tpm) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Nascent RNA") +
  theme(legend.position = "none")


g1.7 <- ggplot(TX.DE.tpm, aes(x = F_log2FC_mTORi, y = F_log2FC_2i, 
                              color = get_dens(F_log2FC_mTORi, F_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TX.DE.tpm) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Total RNA") +
  theme(legend.position = "none")

g1.8 <- ggplot(TX.DE.tpm, aes(x = L_log2FC_mTORi - F_log2FC_mTORi, 
                      y = L_log2FC_2i - F_log2FC_2i,
                      color = get_dens(L_log2FC_mTORi - F_log2FC_mTORi, L_log2FC_2i - F_log2FC_2i, 1000))) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  geom_vline(xintercept = 0, color = "black", lty = 2) +
  geom_point(pch = 15, cex = 0.8) +
  geom_text(x = -3.5, y = 3.5, label = paste("n =", nrow(TX.DE.tpm) ), cex = 4.5, hjust = 0, color = "black") +
  theme_setting +
  scale_color_gradientn(colours = c("#b5b5b5", "#8bc926", "#6fbf89", "#008cf0")) +
  scale_x_continuous(name="log2FC contrast mTORi 1d", breaks=c(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  scale_y_continuous(name="log2FC contrast 2i 2d", breaks=(-3:3), labels=(-3:3), limits = c(-3.5, 3.5)) +
  labs(title = "Nascent RNA - Total RNA") +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1.6, g1.7, g1.8, nrow = 1),
       filename = "Fig2_LRNA_FRNA_log2FC_comparisons.png", path = "figs",
       device = "png", width = 12, height = 4)
