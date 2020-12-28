# Descriptive analysis of TU annotations

# TTseq LRNA bam files processed with TU filter
# with setting:
#   Merge feature: protein_coding, licnRNA
#   STAN method: log poisson
#   Bin: 200 bp
#   Combine: exon

# Rui Shao
# 2020 Feb
# ------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# pacman::p_load(Rsubread)
# run read counting on TU annotation with FeatureCounts
# summarise read distribution and test differential expression
source("F1_src_TU_anno.R")

# TU annotations read abundance
nascent_reads <- read.table("data/Fig1_Tx_read_abundance.txt",
                            header = T, sep = "\t")
nascent_reads <- readRDS("data/nascent_reads.RData")
trans <- function(x){pmin(x, 2) + 0.005*pmax(x-2,0)}
yticks <- c(0, 0.5, 1, 1.5, 2, 100)

ggplot(nascent_reads, aes(x = location, y = trans(Abundance), fill = Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("Nascent RNA Reads Abundance %") +
  scale_fill_manual(values = colors_20[c(13, 13, 13, 2, 2, 20, 20, 7, 7)]) +
  geom_rect(aes(xmin=0, xmax=6, ymin=2, ymax=2.2), fill="white") +
  scale_y_continuous(limits=c(0, NA), breaks=trans(yticks), labels=yticks) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())

ggsave(filename = "Fig1_TU_reads_abundance.pdf", path = "figs",
       device = "pdf", width = 5, height = 4 )

# TPM
readRDS("data/total_TPM.RData") %>%
  ggplot(aes(x = location, y = log10(TPM), fill = Sample)) +
  geom_boxplot(outlier.colour=NA, lwd=0.3) +
  xlab("") +
  ylab("TPM") +
  # scale_y_log10() +
  scale_fill_manual(values = colors_20[c(13, 13, 13, 2, 2, 20, 20, 7, 7)]) +
  scale_y_continuous(limits=c(-1, 3.5), breaks=c(-1, 0, 1, 2, 3), labels=c(0.1, 1, 10, 100, 1000)) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
ggsave(filename = "Fig1_TU_TPM.pdf", path = "figs",
       device = "pdf", width = 5, height = 4 )

# ------------------------------------------------------------------------------
# measure TU interval consistency with jaccard coefficient
jaccard_list <- jaccard_index_pair(TU.list)

g1 <- ggplot(jaccard_list[["protein_coding"]],
       aes(x = Var2, y = Var1, fill = Jaccard)) +
  geom_tile(width=.9, height=.9) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1),
                       values = c(0, .3, .4, .45, .5, .55, .6, .9, 1)) +
  xlab("") + ylab("") + ggtitle("mRNA") +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

g2 <- ggplot(jaccard_list[["uaRNA"]],
       aes(x = Var2, y = Var1, fill = Jaccard)) +
  geom_tile(width=.9, height=.9) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1),
                       values = c(0, .3, .4, .45, .5, .55, .6, .9, 1)) +
  xlab("") + ylab("") + ggtitle("uaRNA") +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

g3 <- ggplot(jaccard_list[["antisense"]],
             aes(x = Var2, y = Var1, fill = Jaccard)) +
  geom_tile(width=.9, height=.9) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1),
                       values = c(0, .3, .4, .45, .5, .55, .6, .9, 1)) +
  xlab("") + ylab("") + ggtitle("asRNA") +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")

g4 <- ggplot(jaccard_list[["intergenic"]],
             aes(x = Var2, y = Var1, fill = Jaccard)) +
  geom_tile(width=.9, height=.9) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 1),
                       values = c(0, .3, .4, .45, .5, .55, .6, .9, 1)) +
  xlab("") + ylab("") + ggtitle("intergenic") +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))

ggsave(plot = grid.arrange(g1, g2, g3, g4, nrow = 1, widths = c(3,3,3,4)),
       filename = "Fig1_TU_location_jaccard_index.png", path = "../figS1/figs",
       device = "png", width = 14, height = 3.5)

# ------------------------------------------------------------------------------
# make MA plots for different classes of TUs
TU.DE.mm9.gr <- readRDS("../fig1/data/TU.DE.mm9.RData")

TU.DE.mm9.gr$location <- as.character(TU.DE.mm9.gr$location)
TU.DE.mm9.gr$location[TU.DE.mm9.gr$location == "protein_coding"] <- "mRNA"
TU.DE.mm9.gr$location[TU.DE.mm9.gr$location == "antisense"] <- "asRNA"

g_list <- list()
for (i in c("mRNA", "intergenic", "uaRNA", "asRNA")) {
  for (x in c("2i", "mTORi")) {
    .title <- paste0(i, " (", sum(TU.DE.mm9.gr$location == i), ")")
    tmp <- with(mcols(TU.DE.mm9.gr)[TU.DE.mm9.gr$location == i, ],
                plot_ma(baseMean = baseMean_LRNA_sp,
                        log2FC = get(paste0("log2FoldChange_LRNA_sp_", x)),
                        title =  ifelse(x == "2i", .title, ""),
                        ylab = "Log2FC",
                        p_col = colors_20[ifelse(x == "2i", 2, 4)]))
    g_list <- c(g_list, list(tmp))
  }
}

ggsave(plot = do.call(grid.arrange, c(g_list, ncol = 2)),
       filename = paste0("Fig1_MAplot_TU_DE_sp_norm.png"), 
       path = "../fig1/figs/",
       device = "png", width = 8, height = 14)

# ------------------------------------------------------------------------------------
# TUs LRNA and FRNA correlation by locations
cor_dat <- data.frame()
for (i in c("mRNA", "intergenic", "uaRNA", "asRNA")) {
  idx <- TU.DE.mm9.gr$location == i
  samples <- c("SL", "2i", "mTORi")
  # samples <- c("SL_1", "SL_2", "2i_1", "2i_2", "mTORi_1", "mTORi_2")
  # na_idx <- tmp.grid[, 1] == tmp.grid[, 2]
  # F_idx <- tmp.grid[, 1] < tmp.grid[, 2]
  
  tmp.grid <- expand.grid(seq_along(samples), seq_along(samples))
  # to_RPK <- function(x) log1p(x) - log(width(TU.DE.mm9.gr)[idx]/1e3)
  to_RPK <- function(x) log1p(rowMeans(x)) - log(width(TU.DE.mm9.gr)[idx]/1e3)
  to_Cor <- function(x) {
    x1 = x[1]; x2 = x[2]
    if (x1 > x2) {
      cor(tmp_LRNA[, x1], tmp_LRNA[, x2])
    } else if (x1 < x2) {
      cor(tmp_FRNA[, x1], tmp_FRNA[, x2])
    } else {
      NA
    }
  }
  
  tmp_LRNA <- cbind(TU.counts.mat.mm9[idx, grep("LRNA.*SL_rep(1|2)", colnames(TU.counts.mat.mm9))] %>% to_RPK(),
                    TU.counts.mat.mm9[idx, grep("LRNA.*2i_2d", colnames(TU.counts.mat.mm9))] %>% to_RPK() ,
                    TU.counts.mat.mm9[idx, grep("LRNA.*mTORi_1d", colnames(TU.counts.mat.mm9))]  %>% to_RPK())
  tmp_FRNA <- cbind(TU.counts.mat.mm9[idx, grep("FRNA.*SL_rep(1|2)", colnames(TU.counts.mat.mm9))] %>% to_RPK(),
                    TU.counts.mat.mm9[idx, grep("FRNA.*2i_2d", colnames(TU.counts.mat.mm9))] %>% to_RPK(),
                    TU.counts.mat.mm9[idx, grep("FRNA.*mTORi_1d", colnames(TU.counts.mat.mm9))] %>% to_RPK())
  
  tmp <- data.frame(x_row = samples[tmp.grid$Var1], 
                    x_col = samples[tmp.grid$Var2],
                    Correlation = apply(tmp.grid, 1, to_Cor))
  tmp$Type <- i
  tmp$RPK <- "Nascent"; tmp$RPK[tmp.grid[, 2] > tmp.grid[, 1]] <- "Total"
  
  cor_dat <- rbind(cor_dat, tmp)
}
cor_dat$x_row <- factor(cor_dat$x_row, levels = rev(samples))
cor_dat$x_col <- factor(cor_dat$x_col, levels = (samples))
cor_dat$Type <- factor(cor_dat$Type, levels = c("mRNA", "intergenic", "uaRNA", "asRNA"))

ggplot(cor_dat, aes(x = x_col, y = x_row)) +
  geom_tile(aes(width=0.9, height=0.9, lty = RPK, fill = Correlation), size = 1) + 
  facet_grid(.~Type, scales = "free", space = "free", switch = "y") +
  xlab("") + ylab("") +
  scale_fill_viridis_c(option = "A", direction = -1, alpha = 0.7, na.value = "grey90") +
  scale_colour_manual(values = c("blue", "green")) +
  scale_size(range = c(0.01, 5), guide = F) +
  theme_setting +
  theme(axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 12, hjust = 1),
        strip.placement = "outside",
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12),
        panel.border = element_blank(), 
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5)) 

ggsave(filename = paste0("Fig1_TU_expr_correlation.png"), 
       path = "../fig1/figs/",
       device = "png", width = 9, height = 3)

# ------------------------------------------------------------------------------------
# explain transcription variation with CpG intensity. scBS-seq
TU.scCpG.mm9 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG/scBS", 
                                    ".bw$", full.names = T),
                         intervals = TU.DE.mm9.gr, fast = T)
colnames(TU.scCpG.mm9) <- gsub(".*(Ser.*|2i.*).cov.bw", "\\1", colnames(TU.scCpG.mm9))
TU.scCpG.mm9 <- TU.scCpG.mm9 / 100

test_CpG_mean_variance <- function(CpG) {
  foreach(i = seq_len(nrow(CpG)), .combine = "rbind") %dopar% {
    x = CpG[i, ]; x = x[!is.na(x)]
    if (length(x) > 0) {
      return(c(mean(x), var(x)))
    } else {
      return(c(NA, NA))
    }
  }
}

registerDoParallel(cores = 10)
scCpG_mat <- cbind(test_CpG_mean_variance(TU.scCpG.mm9[, grep("Ser", colnames(TU.scCpG.mm9))]),
                   test_CpG_mean_variance(TU.scCpG.mm9[, grep("2i", colnames(TU.scCpG.mm9))]))
colnames(scCpG_mat) <- c("mean_sl", "var_sl", "mean_2i", "var_2i")

# bulk bisulfide sequencing
TU.CpG.mm9.gb <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG", 
                                  ".bw$", full.names = T),
                       intervals = TU.DE.mm9.gr,
                       fast = T)[, c(4, 2, 3)]
TU.CpG.mm9.tss <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG", 
                                  ".bw$", full.names = T),
                       intervals = promoters(TU.DE.mm9.gr, upstream = 500, downstream = 0),
                       fast = T)[, c(4, 2, 3)]

# fit first order kinetics to estimate DNA demethylation rate
.time <- cbind(1, c(0, 1, 7))

initial_CpG_and_slope_tss <- solve(t(.time) %*% .time) %*% t(.time) %*% t(log(as.matrix(TU.CpG.mm9.tss))) %>% t()
initial_CpG_and_slope_gb <- solve(t(.time) %*% .time) %*% t(.time) %*% t(log(as.matrix(TU.CpG.mm9.gb))) %>% t()

cls <- cut(TU.CpG.mm9.gb[, 1],
           quantile(TU.CpG.mm9.gb[, 1], seq(0, 1, 0.13), na.rm = T))

# plot DNA methylation change and tx change
dat <- data.frame()
for (i in c("mRNA", "intergenic")) {
  idx <- TU.DE.mm9.gr$location == i
  tmp <- NULL
  for (j in levels(cls)) {
    idx2 = which(idx & (cls == j))
    tmp <- cbind(i, j,
                 c(length(idx2) / sum(idx) * 100, 
                   single_variance_explained(log((TU.DE.mm9.gr$baseMean_LRNA_sp / width(TU.DE.mm9.gr) * 1e3)[idx2]), 
                                             (TU.CpG.mm9.gb - TU.CpG.mm9.tss)[idx2, 2], is.cor = T),
                   single_variance_explained(TU.DE.mm9.gr$log2FoldChange_LRNA_sp_2i[idx2], 
                                             initial_CpG_and_slope_gb[idx2, 2], is.cor = T),
                   median(initial_CpG_and_slope_gb[idx2, 2], na.rm = T)),
                 c("Occurrence(%)", "Cor: Tx RPK ~\nMethyl.(GB - TSS)", "Cor: Tx Log2FC ~\nDemethylation rate", "Demethylation\nrate"))
    dat <- rbind(dat, tmp)
  }
}
colnames(dat) <- c("Type", "Methylation", "Values", "Features")
dat$Values <- as.numeric(as.character(dat$Values))
dat$Methylation <- factor(dat$Methylation, levels = levels(cls))
dat$Features <- factor(dat$Features, 
                       levels = c("Occurrence(%)", "Demethylation\nrate", "Cor: Tx RPK ~\nMethyl.(GB - TSS)", "Cor: Tx Log2FC ~\nDemethylation rate"))

ggplot(dat, aes(x = Methylation, y = Values, color = Type, group = Type)) +
  geom_hline(aes(yintercept=0), color="grey50", size = 0.2) +
  geom_smooth(size = 0, alpha = 0.1, se = T, method = "loess") +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(size = 1) +
  facet_grid(Features~., scales="free_y", switch="y") +
  scale_y_continuous(position = "right") +
  
  scale_color_manual(values = c("tan2", "steelblue")) +
  scale_x_discrete(labels = c(0.1, 0.27, 0.42, 0.54, 0.62, 0.70, 0.78)) +
  ylab("") + xlab("\nTU Methylation level") +
  theme_setting +
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        legend.position = "top")

ggsave(filename = "../fig1/figs/Fig1_TU_methylation_tx.png",
       width = 3.5, height = 6, device = "png")

# ---------------------------------------------------------------------------------




