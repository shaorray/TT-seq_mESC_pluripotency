# Rui Shao 2020 Aug
# DE test
# -----------------------------------------------------------------------------------------------------
# This part includes:
#     1. TU differential expression (LRNA) by sample and TU types
#     2. TU annotation similarity test in FigS2
#     3. LRNA and FRNA correlation by sample and TU types in FigS2
#     4. TU differential expression with evolution conservation
# -----------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# make MA plots for different classes of TUs
TU.DE.mm9.gr <- readRDS("../fig1/data/TU.DE.mm9.RData")

TU.DE.mm9.gr$location <- as.character(TU.DE.mm9.gr$location)
TU.DE.mm9.gr$location[TU.DE.mm9.gr$location == "protein_coding"] <- "mRNA"
TU.DE.mm9.gr$location[TU.DE.mm9.gr$location == "antisense"] <- "asRNA"

bi_other.idx <- findOverlaps(TU.DE.mm9.gr[TU.DE.mm9.gr$location == "intergenic"],
                             bi_intergenic.other.mm9.gr) %>% queryHits() # bi_intergenic.other.mm9.gr from "F1_src_TU_anno.R"

g_list <- list()
for (i in c("mRNA", "intergenic", "uaRNA", "asRNA")) {
  for (x in c("2i", "mTORi")) {
    .title <- paste0(i, " (", sum(TU.DE.mm9.gr$location == i), ")")
    tmp <- with(mcols(TU.DE.mm9.gr)[TU.DE.mm9.gr$location == i, ],
                plot_ma(baseMean = baseMean_LRNA_sp,
                        log2FC = get(paste0("log2FoldChange_LRNA_sp_", x)),
                        title =  ifelse(x == "2i", .title, ""),
                        ylab = "Log2FC",
                        p_col = colors_20[ifelse(x == "2i", 2, 4)],
                        pin_point = .ifelse(i == "intergenic", bi_other.idx, NA) ))
    g_list <- c(g_list, list(tmp))
  }
}

ggsave(plot = do.call(grid.arrange, c(g_list, ncol = 2)),
       filename = paste0("Fig2_MAplot_TU_DE_sp_norm2.png"), 
       path = "../fig2/figs/",
       device = "png", width = 8, height = 14)

# ------------------------------------------------------------------------------
# non coding TU: intergenic differential expression 
readRDS("data/bi_intergenic_DE.RData") %>% # from "F1_src_TU_anno.R"
  ggplot(aes(x = interaction(Enhancer, Sample),
             y = log2FoldChange, fill = Sample)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_boxplot(aes(linetype = Direction), outlier.size = 0, notch=TRUE) +
  ylim(c(-9, 9)) +
  scale_linetype_manual(name = "Direction", values = c(1, 11)) +
  scale_x_discrete(labels = c("2i \nEnhancer", "2i\nOther", "mTORi\nEnhancer", "mTORi\nOther")) +
  scale_fill_manual(values = colors_20[c(2, 20)]) +
  annotate("text", x=1.4, y=9,label = "p_value < 1e-5", cex=4,parse=TRUE) +
  annotate("text", x= c(t(matrix(c(1:4 - .18, 1:4 + .18), nrow = 4))), y=8,
           label= c("ns", "ns", "'*'", "'*'", "ns", "ns", "ns", "ns"), cex=4, parse=TRUE) +
  xlab("") +
  labs(title="Bidirectional Intergenic TU",x="") +
  theme_setting -> g3

readRDS("data/uni_intergenic_DE.RData") %>%
  ggplot(aes(x = interaction(Enhancer, Sample),
             y = log2FoldChange, fill = Sample)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_boxplot(outlier.size = 0, width = 0.5, notch=TRUE) +
  ylim(c(-9, 9)) +
  scale_fill_manual(values = colors_20[c(2, 20)]) +
  scale_x_discrete(labels = c("Enhancer", "Other", "mTORi\nEnhancer", "mTORi\nOther")) +
  annotate("text", x=1.5, y= 9,label = "p_value < 1e-5", cex=4,parse=TRUE) +
  annotate("text", x= 1:4, y= 8, label= rep("ns", 4), cex=4, parse=TRUE) +
  xlab("") +
  labs(title="Unidirectional Intergenic TU",x="") +
  theme_setting +
  theme(legend.position = "none") -> g4

# ChromHMM states overlap coefficient
coef_dat <- readRDS("data/coef_dat.RData")
ggplot(coef_dat, aes(x = chromHMM, y = Estimate, color = V1)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_errorbar(aes(ymin = Estimate - coef_dat$`Std. Error`,
                    ymax = Estimate + coef_dat$`Std. Error`), width = .2,
                position = position_dodge(.5)) +
  geom_point(position = position_dodge(0.5)) +
  scale_color_manual(values = colors_9[c(6, 2)]) +
  coord_flip() +
  labs(color="", title = "log2FC ~ ChromHMM") +
  xlab("") + ylab("coefficient") +
  theme_setting -> g5

ggsave(plot = egg::ggarrange(g4, g3, g5,nrow = 1, widths = c(6,6,4)),
       filename = "Fig2_intergenic_TU_DE.pdf", path = "figs",
       device = "pdf", width = 12, height = 4 )

# ------------------------------------------------------------------------------------
to_RPK <- function(x, idx = NULL) {
  if (is.null(idx)) idx = seq_along(TU.DE.mm9.gr)
  log(rowMeans(x) - width(TU.DE.mm9.gr[idx])/1e3)
}
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

# TUs LRNA correlation by baseMean bin (in log)
zero_expr_idx <- TU.DE.mm9.gr$baseMean_LRNA_sp == 0
TU.DE.mm9.gr$RPK <- (TU.DE.mm9.gr$baseMean_LRNA) / (width(TU.DE.mm9.gr) / 1e3)
rpk_bin <- cut(log(TU.DE.mm9.gr$RPK[!zero_expr_idx]), 
               quantile(log(TU.DE.mm9.gr$RPK[!zero_expr_idx]), seq(0, 1, 0.2)))

TU.counts.mat.mm9_sp_norm <- sweep(TU.counts.mat.mm9[, grep("LRNA", colnames(TU.counts.mat.mm9))],
                                   MARGIN = 2, SizeFactorCal(sp_L_mat[1:3, ]), "/")

TU.RPK.combine_rep.mat.mm9 <- sapply(c("LRNA_SL", "LRNA_2i_2d", "LRNA_mTORi_1d"), 
                                     function(x)
                                       TU.counts.mat.mm9_sp_norm[!zero_expr_idx, 
                                                                 grep(x, colnames(TU.counts.mat.mm9_sp_norm))] %>%  
                                       "/"(width(TU.DE.mm9.gr[!zero_expr_idx]) / 1e3) %>%
                                       rowMeans() %>% log())

rpk_lvl <- quantile(log(TU.DE.mm9.gr$RPK[!zero_expr_idx]), seq(0, 1, 0.1))[(1:5)*2]
names(rpk_lvl) <- levels(rpk_bin)

expr_bin.mat <- NULL
for (i in c("mRNA", "intergenic", "uaRNA", "asRNA")) {
  i.idx <- TU.DE.mm9.gr$location[!zero_expr_idx] == i
  for (j in levels(rpk_bin)) {
    j.idx <- rpk_bin == j
    tmp_dat <- TU.RPK.combine_rep.mat.mm9[j.idx & i.idx, ]
    tmp_dat <- tmp_dat[!apply(tmp_dat, 1, function(x) any(is.infinite(x) | is.na(x))), ]
    expr_bin.mat <- rbind(expr_bin.mat,
                          rbind(c(i, rpk_lvl[j], "2i", "RPK",
                                  cor(tmp_dat[, c(1, 2)], 
                                      use = "pairwise.complete.obs")[1, 2]),
                                c(i, rpk_lvl[j], "mTORi", "RPK",
                                  cor(tmp_dat[, c(1, 3)],
                                      use = "pairwise.complete.obs")[1, 2]))
    )
  }
}
colnames(expr_bin.mat) <- c("Type", "RPK", "Sample", "Expr", "Cor")
expr_bin.mat <- as.data.frame(expr_bin.mat)
expr_bin.mat$RPK <- expr_bin.mat$RPK %>% as.character() %>% as.numeric() %>% round(2)
# expr_bin.mat$RPK <- gsub(".*\\,(.*)]", "\\1", as.character(expr_bin.mat$RPK)) %>% as.numeric() %>% round(1)
# expr_bin.mat$RPK <- quantile(log(TU.DE.mm9.gr$RPK[!zero_expr_idx]), seq(0, 1, 0.1))[(1:5)*2]
expr_bin.mat$Cor <- as.numeric(as.character(expr_bin.mat$Cor))
expr_bin.mat$Type <- factor(expr_bin.mat$Type, levels = c("mRNA", "intergenic", "uaRNA", "asRNA"))

ggplot(expr_bin.mat, aes(x = RPK, y = Cor, color = Type)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_point() +
  geom_line(lwd = 1) +
  # geom_smooth(alpha=0.1, method = "loess", span = 0.95, level=0.90) +
  facet_grid(.~Sample) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  xlab("Mean transcription levels (log RPK)") +
  ylab("Correlation with SL") +
  ylim(c(-0.6, 1)) +
  theme_setting +
  theme(panel.grid.major = element_line(), 
        strip.text = element_text(size = 12),
        legend.position = "top",
        legend.key.size = unit(0.2, "cm") )
ggsave(filename = "Fig2_sample_log_expression_bin_correlation.png",
       path = "../fig2/figs/", width = 5, height = 4, device = "png")


# (Fig EV3) measure TU interval consistency with jaccard coefficient
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
       filename = "FigS2_TU_location_jaccard_index.png", path = "../figS2/figs",
       device = "png", width = 14, height = 3.5)

# ------------------------------------------------------------------------------------
# (Fig EV3) TUs LRNA and FRNA correlation by locations
cor_dat <- data.frame()
for (i in c("mRNA", "intergenic", "uaRNA", "asRNA")) {
  idx <- TU.DE.mm9.gr$location == i
  samples <- c("SL", "2i", "mTORi")
  # samples <- c("SL_1", "SL_2", "2i_1", "2i_2", "mTORi_1", "mTORi_2")
  # na_idx <- tmp.grid[, 1] == tmp.grid[, 2]
  # F_idx <- tmp.grid[, 1] < tmp.grid[, 2]
  
  tmp.grid <- expand.grid(seq_along(samples), seq_along(samples))
  
  tmp_LRNA <- cbind(TU.counts.mat.mm9[idx, grep("LRNA.*SL_rep(1|2)", colnames(TU.counts.mat.mm9))] %>% to_RPK(idx = idx),
                    TU.counts.mat.mm9[idx, grep("LRNA.*2i_2d", colnames(TU.counts.mat.mm9))] %>% to_RPK(idx = idx) ,
                    TU.counts.mat.mm9[idx, grep("LRNA.*mTORi_1d", colnames(TU.counts.mat.mm9))]  %>% to_RPK(idx = idx))
  tmp_FRNA <- cbind(TU.counts.mat.mm9[idx, grep("FRNA.*SL_rep(1|2)", colnames(TU.counts.mat.mm9))] %>% to_RPK(idx = idx),
                    TU.counts.mat.mm9[idx, grep("FRNA.*2i_2d", colnames(TU.counts.mat.mm9))] %>% to_RPK(idx = idx),
                    TU.counts.mat.mm9[idx, grep("FRNA.*mTORi_1d", colnames(TU.counts.mat.mm9))] %>% to_RPK(idx = idx))
  
  tmp <- data.frame(x_row = samples[tmp.grid$Var1], 
                    x_col = samples[tmp.grid$Var2],
                    Correlation = apply(tmp.grid, 1, to_Cor)) # toCor usese global tmp_LRNA and tmp_FRNA in the current loop
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

ggsave(filename = paste0("FigS2_TU_expr_correlation.png"), 
       path = "../figS2/figs/",
       device = "png", width = 9, height = 3)

# ------------------------------------------------------------------------------------
# (Fig EV3) log2FC ~ Evolution conservation + direction

dat_phasrCons <- data.frame(Type = TU.DE.intergenic$enhancer_direction,
                            phastCons = TU.DE.intergenic$mm10.60way.phastCons)
dat_phasrCons$Type <- factor(dat_phasrCons$Type, 
                             c("Unidirectional_Enhancer", "Unidirectional_TX",
                               "Bidirectional_Enhancer", "Bidirectional_TX"))

t.test(dat_phasrCons$phastCons[dat_phasrCons$Type == "Unidirectional_Enhancer"],
       dat_phasrCons$phastCons[dat_phasrCons$Type == "Unidirectional_TX"]) # p-value < 2.2e-16
t.test(dat_phasrCons$phastCons[dat_phasrCons$Type == "Bidirectional_Enhancer"],
       dat_phasrCons$phastCons[dat_phasrCons$Type == "Bidirectional_TX"]) # p-value < 2.2e-16

ggplot(dat_phasrCons, aes(x = Type, y = phastCons, fill = Type), alpha = 0.7) +
  geom_hline(yintercept = median(dat_phasrCons$phastCons[dat_phasrCons$Type == "Unidirectional_Enhancer"]), 
             lty = 2, color = "grey50") +
  geom_violin() +
  geom_boxplot(outlier.color = NA, width = 0.15) +
  ylim(c(-10, 400)) + xlab("") + ylab("phastCons 60way") +
  scale_fill_viridis_d(begin = 0.3, end = 0.95) +
  theme_setting +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "FigS2_phastCons_intergenic_TU_TSS.png", 
       device = "png", path = "../figS2/figs/", width = 4, height = 5)


# ------------------------------------------------------------------------------------
# (not shown) explain transcription variation with CpG intensity. scBS-seq
TU.scCpG.mm9 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG/scBS", 
                                    ".bw$", full.names = T),
                         intervals = TU.DE.mm9.gr, fast = F)
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
                          fast = F)[, c(4, 2, 3)]
TU.CpG.mm9.tss <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG", 
                                      ".bw$", full.names = T),
                           intervals = promoters(TU.DE.mm9.gr, upstream = 500, downstream = 0),
                           fast = F)[, c(4, 2, 3)]

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
                   # single_variance_explained(log((TU.DE.mm9.gr$baseMean_LRNA_sp / width(TU.DE.mm9.gr) * 1e3)[idx2]), 
                   #                           (TU.CpG.mm9.gb - TU.CpG.mm9.tss)[idx2, 2], is.cor = T),
                   single_variance_explained(TU.DE.mm9.gr$log2FoldChange_LRNA_sp_2i[idx2], 
                                             TU.CpG.mm9.gb[idx2, 1] - TU.CpG.mm9.gb[idx2, 3], is.cor = T),
                   single_variance_explained(TU.DE.mm9.gr$log2FoldChange_LRNA_sp_2i[idx2], 
                                             -initial_CpG_and_slope_gb[idx2, 2], is.cor = T),
                   -median(initial_CpG_and_slope_gb[idx2, 2], na.rm = T)),
                 c("Occurrence(%)",
                   # "Cor: Tx RPK ~\nMethyl.(GB - TSS)", 
                   "Cor: Tx Log2FC ~\nDemethylation degree", 
                   "Cor: Tx Log2FC ~\nDemethylation rate", 
                   "Demethylation\nrate"))
    dat <- rbind(dat, tmp)
  }
}
colnames(dat) <- c("Type", "Methylation", "Values", "Features")
dat$Values <- as.numeric(as.character(dat$Values))
dat$Methylation <- factor(dat$Methylation, levels = levels(cls))
dat$Features <- factor(dat$Features, 
                       levels = c("Occurrence(%)",
                                  "Demethylation\nrate",
                                  # "Cor: Tx RPK ~\nMethyl.(GB - TSS)", 
                                  "Cor: Tx Log2FC ~\nDemethylation degree", 
                                  "Cor: Tx Log2FC ~\nDemethylation rate"))

ggplot(dat, aes(x = Methylation, y = Values, color = Type, group = Type)) +
  geom_hline(aes(yintercept=0), color="grey50", size = 0.2) +
  geom_smooth(size = 0, alpha = 0.1, se = T, method = "loess") +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1) +
  facet_grid(Features~., scales="free_y", switch="y") +
  scale_y_continuous(position = "right") +
  
  scale_color_manual(values = c("tan2", "steelblue")) +
  scale_x_discrete(labels = c(0.1, 0.27, 0.42, 0.54, 0.62, 0.70, 0.78)) +
  ylab("") + xlab("\nSerum methylation level") +
  theme_setting +
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 12),
        panel.spacing = unit(1, "lines"),
        legend.position = "top")

ggsave(filename = "figs/Fig2_TU_methylation_tx.png",
       width = 3.5, height = 6, device = "png")

# ---------------------------------------------------------------------------------- #
# (Fig EV3) variance of FRNA log2FC explained by LRNA log2FC
TU.DE.mm9.log2FC <- data.frame(location = TU.DE.mm9.gr$location,
                               LRNA_lfc_2i = TU.DE.mm9.gr$log2FoldChange_LRNA_sp_2i,
                               FRNA_lfc_2i = TU.DE.mm9.gr$log2FoldChange_FRNA_sp_2i,
                               LRNA_lfc_mTORi = TU.DE.mm9.gr$log2FoldChange_LRNA_sp_mTORi,
                               FRNA_lfc_mTORi = TU.DE.mm9.gr$log2FoldChange_FRNA_sp_mTORi) %>%
  dplyr::filter(complete.cases(.))

dat_F_L <- NULL
for (type in c("mRNA", "intergenic", "uaRNA", "conRNA", "asRNA", "daRNA")) {
  tmp <- TU.DE.mm9.log2FC[TU.DE.mm9.log2FC$location == type, ]
  tmp[, 2:5] <- trim_quantile(tmp[, 2:5], 0.95)
  dat_F_L <- rbind(dat_F_L,
                   data.frame(Sample = c("2i 2d", "mTORi 1d"),
                              Type = type,
                              Cor = c(cor(tmp[, 2], tmp[, 3]),
                                      cor(tmp[, 4], tmp[, 5]))
                   ))
}
dat_F_L$Type <- factor(dat_F_L$Type, levels = c("mRNA", "intergenic", "asRNA", "conRNA", "uaRNA", "daRNA"))

ggplot(dat_F_L, aes(x = Type, y = Cor, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
  ylab("log2FC correlation total RNA ~ labeled RNA") + xlab("") +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "../figS2/figs/FigS2_log2FC_correlation_FRNA_LRNA.png",
       width = 6, height = 5, device = "png")

g1 <- TU.DE.mm9.log2FC[TU.DE.mm9.log2FC$location == "mRNA", 2:3] %>%
  ggplot(aes(x = FRNA_lfc_2i, y = LRNA_lfc_2i, 
             color = get_dens(FRNA_lfc_2i, LRNA_lfc_2i, n.grid = 100))) +
  geom_point(cex = 0.5) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.5, vjust = 1.2, 
           label = paste0(" r = ", 0.626, "\nn = ", 12007)) +
  scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.9) +
  xlab("Total RNA log2FC") + ylab("Labeled RNA log2FC") + ggtitle("mRNA (2i_2d vs SL)") +
  xlim(c(-5, 4)) + ylim(c(-5, 4)) +
  theme_setting +
  theme(legend.position = "none")

g2 <- TU.DE.mm9.log2FC[TU.DE.mm9.log2FC$location == "mRNA", 4:5] %>%
  ggplot(aes(x = FRNA_lfc_mTORi, y = LRNA_lfc_mTORi, 
             color = get_dens(FRNA_lfc_mTORi, LRNA_lfc_mTORi, n.grid = 100))) +
  geom_point(cex = 0.5) +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.5, vjust = 1.2, 
           label = paste0(" r = ", 0.64, "\nn = ", 12007)) +
  scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.9) +
  xlab("Total RNA log2FC") + ylab("Labeled RNA log2FC") + ggtitle("mRNA (mTORi_1d vs SL)") +
  xlim(c(-5, 4)) + ylim(c(-5, 4)) +
  theme_setting +
  theme(legend.position = "none")

ggsave(grid.arrange(g1, g2, ncol = 1),
       filename = "../figS2/figs/FigS2_scatter_log2FC_correlation_FRNA_LRNA.png",
       width = 5, height = 7, device = "png")

# ---------------------------------------------------------------------------------- #
# (Fig EV1) mTORi Bulut et al. RNA-seq comparison
Tx_mm9_Bulut <- list.files('/mnt/0E471D453D8EE463/GEO_SL_2i_RNA/output', full.names = T) %>% paste0("/abundance.tsv") %>% 
  SummarizedExperiment::readKallisto(as = 'matrix', what = "est_counts")  %>% 
  as.data.frame() %>%
  dplyr::filter(grepl("^ENS", rownames(.))) %>%
  as.matrix() %>%
  keepOneTx(rowname_gene_id = T, is_gene_sum = T) 

log2FC_mm9_Bulut <- Tx_mm9_Bulut[, 2:1] %>% log2() %>% rowDiffs() %>% `rownames<-`(., rownames(Tx_mm9_Bulut))

# res_FRNA_mTORi is from "F1_src_LoadReadCounts.R"
log2FC_mm9_mTORi_1d <- res_FRNA_mTORi[rank(res_FRNA_mTORi$baseMean) %in% seq(nrow(res_FRNA_mTORi), nrow(res_FRNA_mTORi)-10000) &
                                        rownames(res_FRNA_mTORi) %in% rownames(Tx_mm9_Bulut), ]

log2FC_mm9_Bulut <- log2FC_mm9_Bulut[rownames(log2FC_mm9_mTORi_1d), ]

data.frame(x = log2FC_mm9_mTORi_1d$log2FoldChange, 
           y = log2FC_mm9_Bulut) %>%
  trim_quantile(q = 0.997) %>% `colnames<-`(., c("x", "y")) %>% as.data.frame() %>%
plot_scatter(.xlab = "log2FC FRNA mTORi 1d", .ylab = "log2FC Bulut et al. mTORi 14d",
             xlim = c(-2, 1.5), ylim = c(-2.5, 2.5))
