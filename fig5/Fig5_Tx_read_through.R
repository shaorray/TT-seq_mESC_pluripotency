# Rui Shao 2020 Aug

# -----------------------------------------------------------------------------------------------------
# This part includes:
#     1. calculate termination read-through distance
#     2. estimated speed coverage in 15 kb termination window
#     3. evaluate speed and read-through relationship
#     4. explain read-through with chromatin features
# -----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../util/utils.R')
source('../util/getCoverage.R')

source("F5_src_last_exon_termination_window.R") # get 15 kb termination windows of 10447 genes
source("F5_src_termination_window_coverage.R") # TT-seq and Pol2s5p ChIP bam coverage

# find termination sites using termination window (-5kb, 15kb) expanding the last exons
# include exon coverage 
# load data -----------------------------------------------------------------------------------
terWindow <- importRanges("data/terWindow.mm9.gff3")
gene.last.exon.gr <- importRanges("data/gene.last.exon.mm9.gff3")

TTseq.terWin.cov_norm <- readRDS("data/TTseq.terWin.cov_norm.RData")
Pol2S5p.MINUTE.terWin.cov <- readRDS("data/Pol2S5p.MINUTE.terWin.cov.RData")

# tx termination sites  -----------------------------------------------------------------------------------------------------
if (!file.exists("data/read_thru_table.RData")) {

  get_ter_site <- function(x) {
    if (any(is.na(scale(x)))) return(0)
    which.max(cumsum(scale(log1p(x))))
  }
  
  TTseq.ter.sites <- Reduce(cbind, lapply(TTseq.terWin.cov_norm, 
                                          function(sample.cov) 
                                            unlist(sapply(sample.cov, get_ter_site))))
  colnames(TTseq.ter.sites) <- names(TTseq.terWin.cov_norm)
  
  Pol2S5p.MINUTE.ter.sites <- Reduce(cbind,lapply(Pol2S5p.MINUTE.terWin.cov, 
                                                  function(sample.cov) 
                                                    unlist(sapply(sample.cov, get_ter_site))))
  colnames(Pol2S5p.MINUTE.ter.sites) <- names(Pol2S5p.MINUTE.terWin.cov)
  
  nascent.ter.sites <- Reduce(cbind, lapply(nascent.terWin.cov, 
                                            function(sample.cov) 
                                              unlist(sapply(sample.cov, get_ter_site))))
  colnames(nascent.ter.sites) <- names(nascent.terWin.cov)
  
  pol2.ter.sites <- Reduce(cbind, lapply(pol2.terWin.cov, function(sample.cov) unlist(sapply(sample.cov, get_ter_site))))
  colnames(pol2.ter.sites) <- names(pol2.terWin.cov)
  
  read_thru_table <- cbind(data.frame(gene_id = terWindow$gene_id), # n = 10447 
                           TTseq.ter.sites,
                           Pol2S5p.MINUTE.ter.sites,
                           pol2.ter.sites,
                           nascent.ter.sites)
  saveRDS(read_thru_table, "data/read_thru_table.RData")
}

# plot a read-through example of Hsp90ab1  ---------------------------------------------------
read_thru_table <- readRDS("data/read_thru_table.RData")

.plot_RT_cmp <- function(x, y, xlab, ylab, main = "") {
  ggplot(data = data.frame(x = x, y = y),
         aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(size = 0.3) +
    scale_color_gradientn(colours = rev(colors_n)) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    # theme_minimal() + 
    theme_setting +
    theme(legend.position = "none")
}

g1.1 <- .plot_RT_cmp(read_thru_table$LRNA_2i, read_thru_table$Marta_PROseq,
             xlab = "TT-seq nascent RNA 2i [bp]", ylab = "PRO-seq Lloret et al. 2i [bp]",
             main = "Read-through distance")
g1.2 <- .plot_RT_cmp(read_thru_table$LRNA_2i, read_thru_table$Jesse_PROseq,
              xlab = "TT-seq nascent RNA 2i [bp]", ylab = "PRO-seq Engreitz et al. 2i [bp]")
g1.3 <- .plot_RT_cmp(read_thru_table$LRNA_2i, read_thru_table$P1_Pol5p_2i_CTR_ALL.mm9.fltd.bam,
              xlab = "TT-seq nascent RNA 2i [bp]", ylab = "Pol II-S5p 2i [bp]")
g1.4 <- .plot_RT_cmp(read_thru_table$LRNA_SL, read_thru_table$Flynn_Pol_II_DMSO,
              xlab = "TT-seq nascent RNA serum [bp]", ylab = "GRO-seq Flynn et al. serum [bp]")

ggsave(grid.arrange(g1.1, g1.2, g1.3, g1.4, nrow = 1), 
       filename = "../figS5/figs/FigS5_terSite_methods_comparison.png",
       width = 16, height = 4)

# plot a read-through example of Hsp90ab1  ---------------------------------------------------
id = 9368
ter.data1 = data.frame(Sample = c("SL", "2i", "mTORi"),
                       Site = as.numeric(read_thru_table[id, 2:4]))
g1 <- data.frame(Bp = rep(seq_len(width(terWindow[id])), 3),
           Sample = c(rep("SL", width(terWindow[id])),
                      rep("2i", width(terWindow[id])),
                      rep("mTORi", width(terWindow[id]))),
           Coverage = c(log1p(TTseq.terWin.cov_norm$LRNA_SL[[id]]),
                        log1p(TTseq.terWin.cov_norm$LRNA_2i[[id]]),
                        log1p(TTseq.terWin.cov_norm$LRNA_mTORi[[id]])),
           Average = c(ifelse(log1p(TTseq.terWin.cov_norm$LRNA_SL[[id]]) > mean(log1p(TTseq.terWin.cov_norm$LRNA_SL[[id]])), "Above", "Below"),
                       ifelse(log1p(TTseq.terWin.cov_norm$LRNA_2i[[id]]) > mean(log1p(TTseq.terWin.cov_norm$LRNA_2i[[id]])), "Above", "Below"),
                       ifelse(log1p(TTseq.terWin.cov_norm$LRNA_mTORi[[id]]) > mean(log1p(TTseq.terWin.cov_norm$LRNA_mTORi[[id]])), "Above", "Below"))
           ) %>% 
  mutate(Sample = factor(Sample, levels = c("SL", "2i", "mTORi"))) %>%
  ggplot(aes(x = Bp / 1000, y = Coverage, fill = Average)) +
  geom_bar(stat="identity") +
  facet_grid(Sample~.) + 
  geom_vline(data = ter.data1, aes(xintercept = Site / 1000), lty = 2) +
  scale_fill_manual(values = colors_9[c(8, 3)]) +
  theme_setting +
  ggtitle("Hsp90ab1") + xlab("") + ylab("LRNA coverage") +
  theme(strip.text.y = element_text(size = 14))


g2 <- data.frame(Bp = rep(seq_len(width(terWindow[id])), 3),
           Sample = c(rep("SL", width(terWindow[id])),
                      rep("2i", width(terWindow[id])),
                      rep("mTORi", width(terWindow[id]))),
           Coverage = c(cumsum(scale(log1p(TTseq.terWin.cov_norm$LRNA_SL[[id]]), T, F)),
                        cumsum(scale(log1p(TTseq.terWin.cov_norm$LRNA_2i[[id]]), T, F)),
                        cumsum(scale(log1p(TTseq.terWin.cov_norm$LRNA_mTORi[[id]]), T, F)))
) %>% 
  mutate(Sample = factor(Sample, levels = c("SL", "2i", "mTORi"))) %>%
  ggplot(aes(x = Bp / 1000, y = Coverage / 1000, color = Sample)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = colors_20[(c(13, 2, 7))]) +
  theme_setting +
  xlab("Termination window (Kb)") + ylab("Δ log density") +
  theme(strip.text.y = element_text(size = 14),
        plot.margin = margin(t = 0, r = 0.8, b = 0, l = 0.15, unit = "cm"))

ggsave(plot = grid.arrange(g1, g2, nrow = 2, heights = c(9, 3)),
       filename = "Fig5_terWindow_cov_example.png", 
       path = "../fig5/figs", device = "png", width = 5, height = 7)

g3 <- data.frame(x = (read_thru_table$LRNA_SL + read_thru_table$LRNA_2i ) / 2e3,
           y = (read_thru_table$LRNA_2i - read_thru_table$LRNA_SL) / 1e3
           ) %>%
  ggplot(aes(x = x, y = y)) +
  geom_hex(bins=45) +
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="Mean read through (kb)") +
  scale_y_continuous(name="2i Δread through (kb)") +
  theme(axis.title=element_text(size=16, face="bold"),
        legend.position = "none",
        panel.border = element_blank())

g4 <- data.frame(x = (read_thru_table$LRNA_SL + read_thru_table$LRNA_mTORi ) / 2e3,
                 y = (read_thru_table$LRNA_mTORi - read_thru_table$LRNA_SL) / 1e3
) %>%
  ggplot(aes(x = x, y = y)) +
  geom_hex(bins=45) +
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="Mean read through (kb)") +
  scale_y_continuous(name="mTORi Δread through (kb)") +
  theme(axis.title=element_text(size=16, face="bold"),
        legend.position = "none",
        panel.border = element_blank())

ggsave(plot = grid.arrange(g3, g4, nrow = 2),
       filename = "Fig5_terWindow_SL_2i_mTORi.png", 
       path = "../fig5/figs", device = "png", width = 3.5, height = 7)

# elongation speed ~ read through length  ---------------------------------------------------
elongation_table <- read.table('../fig4/data/elongation_speed_table.txt')

gene_ids <- intersect.Vector(read_thru_table$gene_id, rownames(elongation_table))
mtch <- match(gene_ids, read_thru_table$gene_id)
elongation_table <- elongation_table[gene_ids, ] %>%
  add_column(read_through_SL = read_thru_table$LRNA_SL[mtch],
             read_through_2i = read_thru_table$LRNA_2i[mtch],
             read_through_mTORi = read_thru_table$LRNA_mTORi[mtch]) %>% 
  as.data.frame() %>%
  dplyr::filter(speed > 0.1 & 
                  read_through_SL < 13000 & read_through_2i < 13000 & read_through_mTORi < 13000)

elongation_table$Speed_cls <- cut(log(elongation_table$speed),
                                  breaks = quantile(log(elongation_table$speed), seq(0, 1, 0.2)),
                                  labels = c("Very slow", "Slow", "Medium", "Fast", "Very fast") )

g5 <- ggplot(elongation_table, 
             aes(x = speed,
                 y = read_through_SL / 1000,
                 color = get_dens(log10(speed), read_through_SL / 1000))) +
  geom_point(size = 0.7) +
  # geom_density2d(color = "black", size = 0.4, bins = 10) +
  scale_color_gradientn(colours = rev(colors_n)) +
  scale_x_log10() +
  coord_cartesian(xlim=c(0.1, 7)) +
  annotate(geom = "text", 
           x = 2.6, y = 13,
           hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(with(elongation_table, 
                                                   cor(log10(speed), read_through_SL / 1000)), 3)),
                          "\np < ", formatC(with(elongation_table, 
                                                  cor.test(log10(speed),
                                                           read_through_SL)$p.val), 
                                             format = "e", digits = 0),
                          "\nn = ", nrow(elongation_table))) +
  ylab("Read through (Kb)") +
  xlab("\nMeasured elongation speed (kb/min)\n") +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), legend.position = "none") +
  annotation_logticks(base = 10, sides = "bottom", scaled = T)


dat_speed_change <- readRDS("../fig4/data/dat_speed_change.RData")
mtch <- match(rownames(dat_speed_change), read_thru_table$gene_id)
est_elongation_table <- data.frame(read_through_SL = read_thru_table$LRNA_SL[mtch],
                                   read_through_2i = read_thru_table$LRNA_2i[mtch],
                                   read_through_mTORi = read_thru_table$LRNA_mTORi[mtch],
                                   log10_est_speed = dat_speed_change$v_hat,
                                   log10_speed_change_2i = dat_speed_change$log10FC_2i,
                                   log10_speed_change_mTORi = dat_speed_change$log10FC_mTORi)
est_elongation_table <- est_elongation_table[complete.cases(est_elongation_table) &
                                               est_elongation_table$log10_est_speed > 0 &
                                               est_elongation_table$log10_est_speed < 1.7, ]

g6 <- ggplot(est_elongation_table, 
             aes(x = log10_est_speed,
                 y = read_through_SL / 1000,
                 color = get_dens(log10_est_speed, read_through_SL / 1000))) +
  geom_point(size = 0.7) +
  # geom_density2d(color = "black", size = 0.4, bins = 10) +
  scale_color_gradientn(colours = rev(colors_n)) +
  # scale_x_log10() +
  annotate(geom = "text", 
           x = 1.5, y = 13,
           hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(with(est_elongation_table, cor(log10_est_speed, (read_through_SL))), 3)),
                          "\np < ", formatC(with(elongation_table, 
                                                 cor.test(log10(speed),
                                                          read_through_SL)$p.val), 
                                            format = "e", digits = 0),
                          "\nn = ", nrow(est_elongation_table))) +
  ylab("Read through (Kb)") +
  xlab("\nEst. elongation speed (a.u.)\n") +
  xlim(c(0, 2)) +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), legend.position = "none") +
  annotation_logticks(base = 10, sides = "bottom", scaled = T)

ggsave(plot = grid.arrange(g5, g6, nrow = 2),
       filename = "Fig5_elongation_speed_termination_window.png", path = "figs",
       device = "png", width = 3.5, height = 7)

# elongation speed class ~ ttseq coverage ---------------------------------------------------
names(gene.last.exon.gr) <- gene.last.exon.gr$gene_id
genes.intercect <- intersect.Vector(terWindow[width(terWindow) == 15000]$gene_id,
                                    intersect.Vector(gene.last.exon.gr$gene_id, rownames(elongation_table)))
elongation_table <- elongation_table[genes.intercect, ]

ttseq.cov.exonlast.list = readBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                         pattern = "LRNA.*SL_(rep1|rep2).*bam$", full.names = T),
                                  intervals = gene.last.exon.gr[genes.intercect],
                                  pair_end = T,
                                  stranded = T,
                                  flanks = c(50, 15000),
                                  new_lens = c(50, 50, 200))
ttseq.cov.exonlast = Reduce("+", ttseq.cov.exonlast.list) / 2

mat = aggregate(log1p(ttseq.cov.exonlast[, 51:300]),
                list(elongation_table$Speed_cls),
                median) %>%
  reshape2::melt(id.vars = 'Group.1')

colnames(mat) <- c("Elongation_speed", "Position", "TTseq_coverage")
mat$Position <- gsub("V", "\\2", mat$Position) %>% as.numeric %>% "-"(49)%>% "*"(75) %>% "/"(1000)

c_col <- c("Very slow"=colors_20[11],
           "Slow"=colors_20[12],
           "Medium"=colors_20[13],
           "Fast"=colors_20[14],
           "Very fast"=colors_20[15])

ggplot(mat, aes(x = Position, y = TTseq_coverage + 0.11, fill = Elongation_speed)) +
  geom_bar(stat = "identity", width = 0.076) +
  geom_rect(aes(xmin = -3.7, xmax = 0, ymin = 0, ymax = 0.11),
            fill = "#0000B2", alpha = 0.01, linetype = 0) +
  geom_rect(aes(xmin = 0, xmax = 16, ymin = 0, ymax = 0.11),
            fill = "grey", linetype = 0) +
  geom_rect(aes(xmin = -3.7, xmax = Inf, ymin = 0.11, ymax = 0.15),
            fill = "white", linetype = 0) +
  facet_grid(rows = vars(Elongation_speed)) +
  scale_fill_manual(values = c_col) +
  xlab("Termination window (Kb)") +
  ylab("Nascent RNA coverage") +
  scale_y_continuous(breaks = c(0.15, 1.15), labels = c("0.0", "1.0"),
                     limits = c(0, 1.15), expand = c(0,0)) +
  theme_setting +
  theme(axis.text=element_text(size=10, face = "plain"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.y = element_text(size = 10, face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 0.5, lineend = "butt"),
        axis.ticks.y.left = element_line(size = 0.5),
        panel.spacing = unit(1.5, "lines"))

ggsave(filename = "Fig5_Elongation_speed_class_read_thru_LRNA.png", device = "png",
       width = 3.5, height = 5, path = "../fig5/figs/")

# termination window speed coverage heatmap -----------------------------------------------
speed_cov_SL = speed_cov_2i = speed_cov_mTORi = NULL
fill_zeros = function(x, len = 15000) if (length(x) == len) x else c(x, rep(0, len - length(x))) 

for (i in seq_along(terWindow)) {
  speed_cov_SL = c(speed_cov_SL, 
                   Rle(fill_zeros(TTseq.terWin.cov_norm[[1]][[i]] / (Pol2S5p.MINUTE.terWin.cov[[3]][i] + 1))) )
  speed_cov_2i = c(speed_cov_2i, 
                   Rle(fill_zeros((TTseq.terWin.cov_norm[[2]][[i]] / (Pol2S5p.MINUTE.terWin.cov[[1]][i] + 1)))) )
  speed_cov_mTORi = c(speed_cov_mTORi, 
                      Rle(fill_zeros((TTseq.terWin.cov_norm[[3]][[i]] / (Pol2S5p.MINUTE.terWin.cov[[4]][i] + 1)))) )
}
speed_cov_SL = resizeCov(speed_cov_SL, df = 20, len = 150)
speed_cov_2i = resizeCov(speed_cov_2i, df = 20, len = 150)
speed_cov_mTORi = resizeCov(speed_cov_mTORi, df = 20, len = 150)
saveRDS(list("speed_cov_SL"=speed_cov_SL, "speed_cov_2i"=speed_cov_2i, "speed_cov_mTORi"=speed_cov_mTORi),
        "data/speed_cov_samples.RData")

row_order = order(read_thru_table$LRNA_2i %/% 100, decreasing = F)
row_order = order(rowSums(speed_cov_SL[, 1:20]), decreasing = F)

colors <- colorRampPalette(c("#fcfced", "#bfb49b", "#c76e58", "#b05b3f", "#a62100"))(100)
breaks <- c(seq(0, max(log1p(speed_cov_SL)), length.out = 100), 1e6)

png("figs/Fig5_Speed_coverage_SL.png", width = 1500, height = 500)
image(t(log1p(speed_cov_SL[order(read_thru_table$LRNA_SL %/% 100, decreasing = F), ])),
      xaxt = "none", yaxt = "none", col = colors, breaks = breaks)
axis(side=1, at=c(0, 0.333, 0.667, 1), labels=c(0, 5, 10, 15))
box()
dev.off()
png("figs/Fig5_Speed_terSite_density_SL.png", width = 1500, height = 500)
# plot(density(read_thru_table$LRNA_SL), axes = T, lwd = 15, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
hist(read_thru_table$LRNA_SL, col = "grey40", breaks = 100,
     axes = T, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()

png("figs/Fig5_Speed_coverage_2i.png", width = 1500, height = 500)
image(t(log1p(speed_cov_2i[order(read_thru_table$LRNA_2i %/% 100, decreasing = F), ])),
      xaxt = "none", yaxt = "none", col = colors, breaks = breaks)
axis(side=1, at=c(0, 0.333, 0.667, 1), labels=c(0, 5, 10, 15))
box()
dev.off()
png("figs/Fig5_Speed_terSite_density_2i.png", width = 1500, height = 500)
# plot(density(read_thru_table$LRNA_2i), axes = T, lwd = 15,  main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
hist(read_thru_table$LRNA_2i, col = "grey40", breaks = 100,
     axes = T, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()

png("figs/Fig5_Speed_coverage_mTORi.png", width = 1500, height = 500)
image(t(log1p(speed_cov_mTORi[order(read_thru_table$LRNA_mTORi %/% 100, decreasing = F), ])),
      xaxt = "none", yaxt = "none", col = colors, breaks = breaks)
axis(side=1, at=c(0, 0.333, 0.667, 1), labels=c(0, 5, 10, 15))
box()
dev.off()
png("figs/Fig5_Speed_terSite_density_mTORi.png", width = 1500, height = 500)
# plot(density(read_thru_table$LRNA_mTORi), axes = T, lwd = 15, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
hist(read_thru_table$LRNA_mTORi, col = "grey40", breaks = 100,
     axes = T, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()

png("figs/Fig5_Speed_coverage_scale.png", width = 100, height = 500)
image(x=1, y=seq(0, max(log1p(speed_cov_SL)), length.out=100),
      z=matrix(seq(0, max(log1p(speed_cov_SL)), length.out=100), 1, 100),
      col = colors,
      ylab='', xlab='', xaxt='n', yaxt='n',  axes=T)
axis(4, at=c(0, seq_len(max(log1p(speed_cov_SL)))), 
     labels=c(0, seq_len(max(log1p(speed_cov_SL)))))
box()
dev.off()

# read through speed comparison
dat_speed_cmp <- data.frame(speed = c(colMedians(speed_cov_SL),
                                      colMedians(speed_cov_2i),
                                      colMedians(speed_cov_mTORi)),
                            sample = c(rep("SL", ncol(speed_cov_SL)),
                                       rep("2i", ncol(speed_cov_SL)),
                                       rep("mTORi", ncol(speed_cov_SL)) ),
                            pos = rep(1:150, 3))
dat_speed_cmp$sample <- factor(dat_speed_cmp$sample, c("SL", "2i", "mTORi"))

ggplot(dat_speed_cmp, aes(x = pos, y = speed, color = sample)) +
  geom_line(lwd = 2) +
  scale_x_continuous(name = "Termination window (kb)",
                     breaks = c(0, 50, 100, 150), labels = c(0, 5, 10, 15)) +
  scale_color_manual(values = colors_20[c(13, 2, 7)]) +
  ylab("Est. median speed (a.u.)") +
  annotate(geom = "text", 
           x = 120, y = 5,
           hjust = "left", vjust = "top", label = "n = 10477") +
  theme_setting

ggsave(filename = "Fig5_Termination_window_sample_speed_coverage.png", device = "png",
       width = 5, height = 4, path = "../fig5/figs")

if (T) {
  # read through indexing SL ~ 2i & SL ~ mTORi ------------------------------------------------------
  # longer is defined in range (0.5, 5 Kb)
  # shorter is defined in range (-5, -0.5 Kb)
  # unchanged is defined in range (-0.5, 0.5 Kb)
  
  longer_2i <- with(read_thru_table, LRNA_2i - LRNA_SL > 500 & LRNA_2i - LRNA_SL < 7500)
  shorter_2i <- with(read_thru_table, LRNA_2i - LRNA_SL < -500 & LRNA_2i - LRNA_SL > -7500)
  unchanged_2i <- with(read_thru_table, LRNA_2i - LRNA_SL > -500 & LRNA_2i - LRNA_SL < 500)
  
  longer_mTORi <- with(read_thru_table, LRNA_mTORi - LRNA_SL > 500 & LRNA_mTORi - LRNA_SL < 7500)
  shorter_mTORi <- with(read_thru_table, LRNA_mTORi - LRNA_SL < -500 & LRNA_mTORi - LRNA_SL > -7500)
  unchanged_mTORi <- with(read_thru_table, LRNA_mTORi - LRNA_SL > -500 & LRNA_mTORi - LRNA_SL < 500)
  
  read_thru_changes <- data.frame(row.names = read_thru_table$gene_id)
  read_thru_changes$rt_2i = read_thru_changes$rt_mTORi = NA
  read_thru_changes$rt_2i[longer_2i] = "Longer"
  read_thru_changes$rt_2i[shorter_2i] = "Shorter"
  read_thru_changes$rt_2i[unchanged_2i] = "Unchanged"
  read_thru_changes$rt_mTORi[longer_mTORi] = "Longer"
  read_thru_changes$rt_mTORi[shorter_mTORi] = "Shorter"
  read_thru_changes$rt_mTORi[unchanged_mTORi] = "Unchanged"
  
  png(filename = "figs/Fig5_tx_read_thru_longer.png", width = 500, height = 500)
  plot_Vennerable(list_1 = which(longer_2i), 
                  list_2 = which(longer_mTORi),
                  name_1 = "2i", name_2 = "mTORi", 
                  color_set = c(colors_20[c(2, 6, 20)]), 
                  color_Text = c(colors_20[c(2, 20)]))
  dev.off()
  
  png(filename = "figs/Fig5_tx_read_thru_shorter.png", width = 400, height = 400)
  plot_Vennerable(list_1 = which(shorter_2i), 
                  list_2 = which(shorter_mTORi),
                  name_1 = "2i", name_2 = "mTORi", 
                  color_set = c(colors_20[c(2, 6, 20)]), 
                  color_Text = c(colors_20[c(2, 20)]))
  dev.off()
  
  png(filename = "figs/Fig5_tx_read_thru_unchanged.png", width = 400, height = 400)
  plot_Vennerable(list_1 = which(unchanged_2i), 
                  list_2 = which(unchanged_mTORi),
                  name_1 = "2i", name_2 = "mTORi", 
                  color_set = c(colors_20[c(2, 6, 20)]), 
                  color_Text = c(colors_20[c(2, 20)]))
  dev.off()
}

if (T) {
  # nucleotide frequency on termination sites -----------------------------------------------------------
  names(terWindow) <- terWindow$gene_id
  get_terSite <- function(terSite, TTS = read_thru_table$LRNA_SL) {
    start(terSite) = start(terSite) + ifelse(as.character(strand(terSite)) == "+", 
                                             TTS, width(terWindow) - TTS)
    end(terSite) = start(terSite) + 1
    terSite
  }
  
  terSite_SL <- get_terSite(terSite = terWindow[read_thru_table$gene_id], TTS = read_thru_table$LRNA_SL)
  terSite_2i <- get_terSite(terSite = terWindow[read_thru_table$gene_id], TTS = read_thru_table$LRNA_2i)
  terSite_mTORi <- get_terSite(terSite = terWindow[read_thru_table$gene_id], TTS = read_thru_table$LRNA_mTORi)
  
  get_nt_freq <- function(intervals) { 
    A_coverage <- interval_pattern_hits(intervals = intervals, pattern = "A", which_genome = "mm9", 
                                        to_coverage = T, out_width = 200)
    T_coverage <- interval_pattern_hits(intervals = intervals, pattern = "T", which_genome = "mm9", 
                                        to_coverage = T, out_width = 200)
    C_coverage <- interval_pattern_hits(intervals = intervals, pattern = "C", which_genome = "mm9", 
                                        to_coverage = T, out_width = 200)
    G_coverage <- interval_pattern_hits(intervals = intervals, pattern = "G", which_genome = "mm9", 
                                        to_coverage = T, out_width = 200)
    list("A_coverage" = A_coverage, "T_coverage" = T_coverage, 
         "C_coverage" = C_coverage, "G_coverage" = G_coverage)
  }
  nt_frequency_SL.list <- get_nt_freq(intervals = terSite_SL + 99)
  nt_frequency_2i.list <- get_nt_freq(intervals = terSite_2i + 99)
  nt_frequency_mTORi.list <- get_nt_freq(intervals = terSite_mTORi + 99)
  
  dat_GC_AT <- cbind(rep(c("SL", "2i", "mTORi"), each = 600),
                     rbind(with(nt_frequency_SL.list, 
                                aggregate(G_coverage + C_coverage - A_coverage - T_coverage,
                                          list(read_thru_changes$rt_mTORi),
                                          mean) %>%
                                  reshape2::melt(id.vars = 'Group.1')),
                           
                           with(nt_frequency_2i.list, 
                                aggregate(G_coverage + C_coverage - A_coverage - T_coverage,
                                          list(read_thru_changes$rt_2i),
                                          mean) %>%
                                  reshape2::melt(id.vars = 'Group.1')),
                           
                           with(nt_frequency_mTORi.list, 
                                aggregate(G_coverage + C_coverage - A_coverage - T_coverage,
                                          list(read_thru_changes$rt_mTORi),
                                          mean) %>%
                                  reshape2::melt(id.vars = 'Group.1')) )
  )
  colnames(dat_GC_AT) <- c("Sample", "Read through", "Position", "Freq")
  dat_GC_AT$Position <- gsub("V", "\\2", dat_GC_AT$Position) %>% as.numeric() %>% "-"(100)
  dat_GC_AT$Sample <- factor(dat_GC_AT$Sample, levels = c("SL", "2i", "mTORi"))
  
  ggplot(dat_GC_AT, aes(x = Position, y = Freq, color = Sample)) +
    geom_line(lwd = 1) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", lwd = 0.5) +
    # geom_point(cex = 0.5) +
    # geom_smooth(method = "loess", span = 0.2, method.args = list(degree=1)) + 
    facet_grid(.~`Read through`) + 
    scale_color_manual(values = colors_20[c(13, 2, 7)]) +
    theme_setting +
    xlab("Position from termination sites (nt)") +
    ylab("(G/C) - (A/T)") +
    theme(strip.text.x = element_text(size = 12, face = "bold"), 
          legend.position = "none")
  ggsave(filename = "Fig5_Termination_site_nt_frequency.png", device = "png",
         width = 11, height = 2.5, path = "figs")
}


if (T) {
  # estimated elongation speed changes ~ ttseq coverage and read through changes ---------------------------------------------------
  genes.intercect <- intersect.Vector(names(gene.last.exon.gr), rownames(read_thru_changes))
  read_thru_changes <- read_thru_changes[genes.intercect, ]
  # read coverages
  ttseq.cov.exonlast.SL = readBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                         pattern = "LRNA.*SL_(rep1|rep2).*bam$", full.names = T),
                                  intervals = gene.last.exon.gr[genes.intercect],
                                  pair_end = T,
                                  stranded = T,
                                  flanks = c(50, 15000),
                                  new_lens = c(50, 50, 200))
  ttseq.cov.exonlast.SL = Reduce("+", ttseq.cov.exonlast.SL) / length(ttseq.cov.exonlast.SL)
  
  ttseq.cov.exonlast.2i = readBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                         pattern = "LRNA_2i_2.*bam$", full.names = T),
                                  intervals = gene.last.exon.gr[genes.intercect],
                                  pair_end = T,
                                  stranded = T,
                                  flanks = c(50, 15000),
                                  new_lens = c(50, 50, 200))
  ttseq.cov.exonlast.2i = Reduce("+", ttseq.cov.exonlast.2i) / length(ttseq.cov.exonlast.2i)
  
  ttseq.cov.exonlast.mTORi = readBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                            pattern = "LRNA.*mTORi_1.*bam$", full.names = T),
                                     intervals = gene.last.exon.gr[genes.intercect],
                                     pair_end = T,
                                     stranded = T,
                                     flanks = c(50, 15000),
                                     new_lens = c(50, 50, 200))
  ttseq.cov.exonlast.mTORi = Reduce("+", ttseq.cov.exonlast.mTORi) / length(ttseq.cov.exonlast.mTORi)
  
  # normalise coverage by exon size factors
  
  sf <- SizeFactorCal(cbind(rowSums(ttseq.cov.exonlast.SL[, 51:100]), # normalize with exon coverages
                            rowSums(ttseq.cov.exonlast.2i[, 51:100]),
                            rowSums(ttseq.cov.exonlast.mTORi[, 51:100])))
  
  dat_rt_cls_cov <- cbind(Sample = rep(c("SL", "2i", "mTORi"), each = 750),
                          rbind(aggregate(log1p(ttseq.cov.exonlast.SL[, 51:300] / sf[1]),
                                          list(read_thru_changes$rt_mTORi),
                                          median) %>%
                                  reshape2::melt(id.vars = 'Group.1'),
                                aggregate(log1p(ttseq.cov.exonlast.2i[, 51:300] / sf[2]),
                                          list(read_thru_changes$rt_2i),
                                          median) %>%
                                  reshape2::melt(id.vars = 'Group.1'),
                                aggregate(log1p(ttseq.cov.exonlast.mTORi[, 51:300] / sf[3]),
                                          list(read_thru_changes$rt_mTORi),
                                          median) %>%
                                  reshape2::melt(id.vars = 'Group.1'))
  )
  
  colnames(dat_rt_cls_cov) <- c("Sample", "Read_thru_changes", "Position", "TTseq_coverage")
  dat_rt_cls_cov$Position <- gsub("V", "\\2", dat_rt_cls_cov$Position) %>% as.numeric %>% "-"(49)%>% "*"(75) %>% "/"(1000)
  
  c_col <- c("SL"=colors_20[13],
             "2i"=colors_20[2],
             "mTORi"=colors_20[7])
  
  ggplot(dat_rt_cls_cov, aes(x = Position, y = TTseq_coverage + 0.15, color = Sample)) +
    geom_line(lwd = 1.2) +
    geom_rect(aes(xmin = -3.7, xmax = 0, ymin = 0, ymax = 0.11),
              fill = "#0000B2", alpha = 0.01, linetype = 0) +
    geom_rect(aes(xmin = 0, xmax = 16, ymin = 0, ymax = 0.11),
              fill = "grey", linetype = 0) +
    geom_rect(aes(xmin = -3.7, xmax = Inf, ymin = 0.11, ymax = 0.15),
              fill = "white", linetype = 0) +
    facet_grid(cols = vars(Read_thru_changes)) +
    scale_color_manual(values = c_col) +
    xlab("Termination window (Kb)") +
    ylab("LRNA coverage") +
    scale_y_continuous(breaks = c(0.15, 1.15), labels = c("0.0", "1.0"),
                       limits = c(0, 1.3), expand = c(0,0)) +
    theme_setting +
    theme(axis.text=element_text(size=12, face = "plain"),
          axis.title=element_text(size=14,face="bold"),
          strip.text.y = element_text(size = 12, face="bold"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.y.left = element_line(size = 0.5, lineend = "butt"),
          axis.ticks.y.left = element_line(size = 0.5),
          panel.spacing = unit(1.5, "lines"),
          strip.text.x = element_text(size = 12, face = "bold"))
  ggsave(filename = "Fig5_Read_thru_changes_coverage.png", device = "png",
         width = 11, height = 2.5, path = "figs")
}

if (T) {
  # read-through changes ~ speed changes -----------------------------------------------------------------------
  genes.idx <- intersect.Vector(rownames(read_thru_table), rownames(estimated_speed))
  
  dat_speed_RT <- cbind(10^(estimated_speed[genes.idx, 2:3]) - 10^(estimated_speed[genes.idx, 1]), 
                        (read_thru_speed_2i - read_thru_speed_SL)[genes.idx],
                        (read_thru_speed_mTORi - read_thru_speed_SL)[genes.idx],
                        read_thru_table[genes.idx, 3:4] - read_thru_table[genes.idx, 2])
  
  colnames(dat_speed_RT) <- c("delta_GB_speed_2i", "delta_GB_speed_mTORi",
                              "delta_RT_speed_2i", "delta_RT_speed_mTORi",
                              "delta_RT_2i", "delta_RT_mTORi")
  
  # gene body speed changes 2i
  idx_2i <- with(dat_speed_RT,
                 delta_RT_speed_2i > (-20) & delta_RT_speed_2i < 10 &
                   # delta_GB_speed_2i > (-5) & delta_GB_speed_2i < 3 &
                   abs(delta_RT_2i) < 6000 & abs(delta_RT_2i) > 500) & 
    complete.cases(dat_speed_RT) & abs(rowSums(dat_speed_RT)) < Inf
  
  cor_2i <- cor.test(dat_speed_RT$delta_RT_speed_2i[idx_2i], 
                     dat_speed_RT$delta_RT_2i[idx_2i])
  
  g3.1 <- ggplot(dat_speed_RT[idx_2i, ], 
                 aes(x = delta_RT_speed_2i, y = delta_RT_2i, 
                     color = get_dens(delta_RT_speed_2i, delta_RT_2i))) +
    geom_point(size = 0.5) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -500, ymax = 500),
              fill = "grey90", alpha = 0.01, linetype = 0) +
    annotate(geom = "text", x = -20, y = 0, hjust = "left", 
             label = paste0("Unchanged ", round(sum(with(dat_speed_RT, abs(delta_RT_2i) < 6000 & abs(delta_RT_2i) > 500)) / nrow(dat_speed_RT) * 100), "%")) +
    annotate(geom = "text",
             x = -20, y = 6000,
             hjust = "left", vjust = "top",
             label = paste(paste("r =", round(cor_2i$estimate, 3)),
                           "\np <", formatC(cor_2i$p.value, format = "e", digits = 2),
                           "\nn =", sum(idx_2i))) +
    scale_color_gradientn(colours = add.alpha(rev(colors_n), 0.5)) +
    xlab("ΔRead-through speed (a.u.)") + ylab("ΔRead-through (bp)") + ggtitle("2i changes") +
    theme_setting +
    theme(legend.position = "none")
  
  # gene body speed changes mTORi
  idx_mTORi <- with(dat_speed_RT,
                    delta_RT_speed_mTORi > (-20) & delta_RT_speed_mTORi < 10 &
                      abs(delta_RT_mTORi) < 6000 & abs(delta_RT_mTORi) > 500) & 
    complete.cases(dat_speed_RT) & abs(rowSums(dat_speed_RT)) < Inf
  
  cor_mTORi <- cor.test(dat_speed_RT$delta_RT_speed_mTORi[idx_mTORi], 
                        dat_speed_RT$delta_RT_mTORi[idx_mTORi])
  
  g3.2 <- ggplot(dat_speed_RT[idx_mTORi, ], 
                 aes(x = delta_RT_speed_mTORi, y = delta_RT_mTORi, 
                     color = get_dens(delta_RT_speed_mTORi, delta_RT_mTORi))) +
    geom_point(size = 0.5) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -500, ymax = 500),
              fill = "grey90", alpha = 0.01, linetype = 0) +
    annotate(geom = "text", x = -20, y = 0, hjust = "left", 
             label = paste0("Unchanged ", round(sum(with(dat_speed_RT, abs(delta_RT_mTORi) < 6000 & abs(delta_RT_mTORi) > 500)) / nrow(dat_speed_RT) * 100), "%")) +
    annotate(geom = "text",
             x = -20, y = 6000,
             hjust = "left", vjust = "top",
             label = paste(paste("r =", round(cor_mTORi$estimate, 3)),
                           "\np <", formatC(cor_mTORi$p.value, format = "e", digits = 2),
                           "\nn =", sum(idx_mTORi))) +
    scale_color_gradientn(colours = add.alpha(rev(colors_n), 0.5)) +
    xlab("\nΔRead-through speed (a.u.)") + ylab("ΔRead-through (bp)") +ggtitle("mTORi changes") +
    theme_setting +
    theme(legend.position = "none")
  
  
  ggsave(plot = grid.arrange(g3.1, g3.2, nrow = 2),
         filename = paste0("FigS6_read_through_gene_body_speed_change.png"), 
         path = "../figS6/figs/",
         device = "png", width = 4, height = 7)
}

# read through explanation with estimated elongation speed, gene length, tx length, etc.-----------------
estimated_speed <- readRDS("../fig4/data/est_speed_mat.RData")
read_thru_table <- readRDS("data/read_thru_table.RData")

rownames(read_thru_table) <- read_thru_table$gene_id
# genes.intercect <- intersect.Vector(rownames(read_thru_table), rownames(estimated_speed))

terWindow_SL <- promoters(terWindow, upstream = 0, downstream = read_thru_table$LRNA_SL + 1)

read_thru_GC = (interval_pattern_hits(intervals = terWindow_SL, pattern = "C", which_genome = "mm9", 
                      to_coverage = F) +
                  interval_pattern_hits(intervals = terWindow_SL, pattern = "G", which_genome = "mm9", 
                                        to_coverage = F)
                  ) / width(terWindow_SL)
names(read_thru_GC) <- read_thru_table$gene_id

sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$Sample == "SL", ]
rownames(sample_Tx_counts) <- sample_Tx_counts$gene_id
sample_Tx_counts$synthesis <- sample_Tx_counts[, c("Copy", "Labeled_rate")] %>% log10() %>% rowSums() %>% "-"(log10(5))

# two_mers <- paste0(c("A", "T", "G", "C")[rep(1:4, 4)], c("A", "T", "G", "C")[rep(1:4, each = 4)])
# a2 = sapply(two_mers, interval_pattern_hits, 
#             intervals = terWindow, weight_pos = read_thru_table$LRNA_SL, which_genome = "mm9")
# a2 = a2 / width(terWindow[rownames(read_thru_table)])
# rownames(a2) <- terWindow$gene_id

# extract speed falls --------------------------------------------------------------------------------
# get_slopes <- function(cov, site) coef(lm(cov[seq_len(site)] ~ seq_len(site) ))[2]
# get_time <- function(cov, site) sum(1/cov[seq_len(site)])
# slopes_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
#                     function(x) get_slopes(cov = log1p(speed_cov_SL[x, ]), 
#                                            site = read_thru_table$LRNA_SL[x] %/% 100 + 10))
# names(slopes_SL) <- read_thru_table$gene_id
# 
# time_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
#                   function(x) get_time(cov = speed_cov_SL[x, ], 
#                                        site = read_thru_table$LRNA_SL[x] %/% 100 + 10))
# names(time_SL) <- read_thru_table$gene_id

# get read through speed --------------------------------------------------------------------------------
read_thru_speed_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
                             function(x) 
                               mean(speed_cov_SL[x, seq_len(ceiling(read_thru_table$LRNA_SL[x] / 100))]))

read_thru_speed_2i <- sapply(seq_len(nrow(speed_cov_2i)), 
                             function(x) 
                               mean(speed_cov_2i[x, seq_len(ceiling(read_thru_table$LRNA_2i[x] / 100))]))

read_thru_speed_mTORi <- sapply(seq_len(nrow(speed_cov_mTORi)), 
                                function(x) 
                                  mean(speed_cov_mTORi[x, seq_len(ceiling(read_thru_table$LRNA_mTORi[x] / 100))]))

names(read_thru_speed_SL) <- names(read_thru_speed_2i) <- names(read_thru_speed_mTORi) <- read_thru_table$gene_id


# GC content from ensembl reference database
library(biomaRt)
gene.list <- getBM(attributes = c("ensembl_gene_id", "transcript_length", "percentage_gene_gc_content"), 
                  mart = useMart("mmusculus_gene_ensembl", biomart = "ensembl"))
gene.list <- gene.list[!duplicated(gene.list$ensembl_gene_id), ]
rownames(gene.list) <- gene.list$ensembl_gene_id

genes.intercect <- intersect.Vector(rownames(estimated_speed), gene.gr$gene_id)

dat_test <- data.frame("Read through length" = read_thru_table[genes.intercect, 2],
                       "Synthesis rate" = sample_Tx_counts[genes.intercect, "synthesis"], # log scale
                       "Gene length" = log10(width(gene.gr[genes.intercect])), # gene.gr is from "F4_src_Tx_elongation.R"
                       "Transcript length" = log10(gene.list[genes.intercect, ]$transcript_length),
                       "G.B. GC%" = gene.list[genes.intercect, ]$percentage_gene_gc_content,
                       "R.T. GC%" = read_thru_GC[genes.intercect] * 100,
                       "Est. G.B. speed" = estimated_speed[genes.intercect, 1], # log10 scale
                       "Est. R.T. speed" = log1p(read_thru_speed_SL[genes.intercect]) )
rownames(dat_test) <- genes.intercect
dat_test <- dat_test[complete.cases(dat_test) & dat_test[, 1] < 15000 &
                       !is.infinite(dat_test$Synthesis.rate) & !is.infinite(dat_test$Est..G.B..speed), ]

get_z_score <- function(mat) scale(log1p(mat))

# chromatin remodeler
mat_ter_cdh <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou/", 
                                                pattern = ".*bw", full.names = T), 
                          intervals = terWindow_SL, fast = T) %>% trim_quantile()
sample_names <- gsub(".*/", "\\2", list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou/", 
                                              pattern = ".*bw", full.names = T))
colnames(mat_ter_cdh) <- sample_names
mat_ter_cdh <- sapply(unique(gsub("_replicate.*", "", sample_names)),
                      function(x) {
                        idx <- grep(x, sample_names)
                        if (length(idx) > 1) {
                          rowMeans(mat_ter_cdh[, idx])
                        } else {
                          mat_ter_cdh[, idx]
                        }
                      }) %>% get_z_score()
rownames(mat_ter_cdh) <- terWindow_SL$gene_id

# DNA methylation
mat_ter_CpG <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG/", 
                                              pattern = ".*bw", full.names = T), 
                        intervals = terWindow_SL, fast = F) %>% trim_quantile()
colnames(mat_ter_CpG) <- gsub(".*DNA_CpG_E14_(.*).bw", "\\1", list.files("/mnt/0E471D453D8EE463/GEO_bw/CpG/", pattern = ".*bw"))
mat_ter_CpG <- get_z_score(mat_ter_CpG)
rownames(mat_ter_CpG) <- terWindow_SL$gene_id

# Accessibility
mat_ter_acc <- cbind(.countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/DHS/", 
                                              pattern = ".*", full.names = T), 
                        intervals = terWindow_SL, fast = F) %>% rowMeans(),
                     .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/FAIRE", 
                                                    pattern = ".*", full.names = T), 
                              intervals = terWindow_SL, fast = F) %>% rowMeans()) %>% 
              trim_quantile()
colnames(mat_ter_acc) <- c("DHS", "FAIRE")
mat_ter_acc <- get_z_score(mat_ter_acc)
rownames(mat_ter_acc) <- terWindow_SL$gene_id

# Domain
mat_ter_domain <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/domain/", 
                                              pattern = "(LaminB|HP1).*bw", full.names = T), 
                        intervals = terWindow_SL, fast = F) %>% trim_quantile()

sample_names <- gsub(".*/", "\\2", list.files("/mnt/0E471D453D8EE463/GEO_bw/domain", 
                                              pattern = "(LaminB|HP1).*bw", full.names = T))
colnames(mat_ter_domain) <- sample_names
mat_ter_domain <- sapply(unique(gsub(".*(Ostapcuk_|ESC_)(.*)(_wt|_ChIP).*", "\\2", sample_names)),
                      function(x) {
                        idx <- grep(x, sample_names)
                        if (length(idx) > 1) {
                          rowMeans(mat_ter_domain[, idx])
                        } else {
                          mat_ter_domain[, idx]
                        }
                      }) %>% get_z_score()

rownames(mat_ter_domain) <- terWindow_SL$gene_id

# RNA associated levels
mat_ter_m6A <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/m6A/", 
                                                 pattern = "SySy.*bw", full.names = T), 
                           intervals = terWindow_SL, fast = F) %>% trim_quantile()
mat_ter_m6A <- sweep(mat_ter_m6A, 2, SizeFactorCal(mat_ter_m6A), "/")
mat_ter_m6A <- log2(rowMeans(mat_ter_m6A[, 1:2]) / rowMeans(mat_ter_m6A[, 3:4])) %>% 
               trim_quantile() %>% scale()
colnames(mat_ter_m6A) <- "m6A"
rownames(mat_ter_m6A) <- terWindow_SL$gene_id

mat_ter_KAS <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/ssDNA/", 
                                              pattern = ".*bw", full.names = T), 
                        intervals = terWindow_SL, fast = F) %>% 
               trim_quantile() %>% rowMeans() %>% get_z_score()
colnames(mat_ter_KAS) <- "KAS"
rownames(mat_ter_KAS) <- terWindow_SL$gene_id

# histone modifications
mat_ter_histone <- .countBW(bw_files = list.files("/mnt/0E471D453D8EE463/GEO_bw/Histone", 
                                                  pattern = ".*bw", full.names = T), 
                            intervals = terWindow_SL, fast = F) %>% trim_quantile()
colnames(mat_ter_histone) <- gsub("_ChIP.*|_rep.*|_wt.*|_WT", "\\1",
                                  gsub(".*(20|GSM).*_(H.*).(mm9|WT.rep1).bw", "\\2", 
                                       list.files("/mnt/0E471D453D8EE463/GEO_bw/Histone", 
                                                  pattern = ".*bw")))

mat_ter_histone <- data.frame(mat_ter_histone) %>% get_z_score()

# combine features
gene.idx <- match(rownames(dat_test), terWindow_SL$gene_id)
dat_test[, -1] <- scale(dat_test[, -1])
dat_test <- cbind(dat_test, 
                  mat_ter_cdh[gene.idx, ],
                  mat_ter_acc[gene.idx, ],
                  mat_ter_domain[gene.idx, ],
                  mat_ter_m6A[gene.idx, ],
                  mat_ter_KAS[gene.idx, ],
                  mat_ter_CpG[gene.idx, 4],
                  mat_ter_histone[gene.idx, ])

colnames(dat_test)[c(5:18, 25:27)] <- c("GB_GC", "RT_GC", "Est.GB_speed", "Est.RT_speed",
                                        "Smarca4", "Chd1", "Chd2", "Chd4", "Chd6", "Chd8", 
                                        "Chd9", "ATAC", "Control", "MNase",
                                        "m6A", "KAS", "RT_DNAme")
# saveRDS(dat_test, "../fig5/data/read_through_dat_test.RData")

# linear explanation
multi_variance_explained(dat_test[,-1], dat_test[,1]) %>% sum() # 0.3058384

dat_test_linear <- dat_test[dat_test$Read.through.length > 100,
                            !grepl("Control|MNase|H3K18ac|H3K18cr|H2BK20ac",
                                  colnames(dat_test)) & !duplicated(colnames(dat_test))]

# view features with non-negative R2
dat_test_linear <- dat_test_linear[, c("Read.through.length", "GB_GC", "RT_GC", 
                                       "Chd4", "Chd6", "Chd1", "Chd8",
                                       "H3.3", "H2AUb", "H3K79me2", "H2AZ", 
                                       "H3K27me3", "RT_DNAme", "LaminB", 
                                       "H3K36me3", "Gene.length", "Est.GB_speed", "Est.RT_speed")]

dat_test_linear$RT_class <- cut(dat_test_linear$Read.through.length,
                         breaks = c(0, 3000, 5000, 7500, 15000), 
                         labels = paste("<", c(3, 5, 7.5, 15)))

dat_test_linear$RT_class <- cut(dat_test_linear$Read.through.length,
                         breaks = c(0, 2000, 3000, 4000, 5000, 6000, 9000, 15000), 
                         labels = paste("<", c(2, 3, 4, 5, 6, 9, 15)))

dat_RT_bin <- data.frame()
for(x in levels(dat_test_linear$RT_class)) {
  dat_RT_bin <- rbind(dat_RT_bin,
                      data.frame(RT_class = x, 
                                 z_score = colMeans(as.matrix(dat_test_linear[dat_test_linear$RT_class == x, -c(1, 19)])), 
                                 feature = colnames(dat_test_linear)[-c(1, 19)]) )
} 
dat_RT_bin$feature <- factor(dat_RT_bin$feature, levels = colnames(dat_test_linear)[-c(1, 19)])

ggplot(dat_RT_bin, aes(x = RT_class, y = feature)) +
  geom_tile(aes(fill = z_score)) +
  scale_fill_viridis_c(name = "Z scores", option = "A", begin = 0.1, end = 0.9) +
  xlab("Read through distance [kb]") + ylab("") + ggtitle("Read through window coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.ticks = element_line(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12)) 

ggsave(filename = "Fig5_Read_through_distance_feature_z_score.png",
       path = "../fig5/figs", device = "png", width = 5, height = 7) 

# linear correlation and R-squared
dat_variance_explained <- data.frame(Feature = rep(colnames(dat_test_linear)[-c(1, 19)], 2),
                                     Value = c(multi_variance_explained(dat_test_linear[, -c(1, 19)], 
                                                                            dat_test_linear$Read.through.length, 
                                                                            is.cor = T),
                                                     multi_variance_explained(dat_test_linear[, -c(1, 19)], 
                                                                              dat_test_linear$Read.through.length, 
                                                                              is.cor = F)),
                                     Type = rep(c("Correlation", "R-squared"), each = 17))
dat_variance_explained$Feature <- factor(dat_variance_explained$Feature, levels = colnames(dat_test_linear)[-c(19)])
sum(dat_variance_explained$Value[dat_variance_explained$Type == "R-squared"]) # 0.2533827 variance explained

ggplot(dat_variance_explained, aes(x = Feature, y = Value, fill = Feature)) +
  geom_bar(stat="identity") +
  facet_grid(Type~., scales="free_y") +
  xlab("") + ylab("") +
  geom_hline(yintercept = 0, color = "grey50") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, hjust = 1),
        axis.ticks = element_line(), 
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12))

ggsave(filename = "FigS5_Read_through_distance_feature_explained.png",
       path = "../figS5/figs", device = "png", width = 5, height = 5) 

# non-linear explanation of read through distance ------------------------------------------------------------------
dat_test_gbm <- dat_test[dat_test$Read.through.length > 100,
                         !grepl("Control|MNase|H3K18ac|H3K18cr|H2BK20ac",
                            colnames(dat_test)) & !duplicated(colnames(dat_test))]

source("../fig5/F5_src_read_through_gbm.R")

# termination sites epigenome feature coverage -----------------------------------------------------------------------
library(abind)

combine_rep_array <- function(cov_array) {
  out_array <- cov_array[, , 1]
  for (i in seq_len(dim(cov_array)[3])[-1]) {
    out_array <- out_array + cov_array[, , i]
  }
  out_array / dim(cov_array)[3]
}

# accessibility
cov_FAIRE <- convertCoverage(list.files("/mnt/0E471D453D8EE463/GEO_bw/FAIRE", "*", full.names = T),
                             intervals = terSite_SL + 1159, bin_width = 10, new_len = 29, 
                             sample_names = rep("FAIRE", 2)) %>% 
              combine_rep_array()

# domain
cov_HP1a <- convertCoverage(list.files("/mnt/0E471D453D8EE463/GEO_bw/domain", "HP1a", full.names = T),
                            intervals = terSite_SL + 1159, bin_width = 10, new_len = 29) %>% 
            combine_rep_array()

# TTseq
cov_TTseq <- .coverBam(list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                     pattern = "LRNA.*(SL_).*bam$", full.names = T),
          intervals = terSite_SL + 1159, paired.end = "extend", df = 10, len = 29) %>%
          Reduce(f = "+", x = .)

# Pol2
cov_Pol2S5p <- .coverBam(list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                    pattern = "P1.*5p.*SL_CTR_ALL.*bam$", full.names = T),
                         intervals = terSite_SL + 1159, paired.end = "ignore", df = 10, len = 29) 


plot_coverage_profile <- function(cov_mat, line.color = "grey50", 
                                  x_left_lab = "-1", x_mid_lab = "0", x_right_lab = "1",
                                  xlab = "", ylab = "", .title = "") {

  cov_mat <- apply(cov_mat[complete.cases(cov_mat), ], 1,
                   function(x) smooth.spline(x = seq_along(x), y = as.numeric(x), df = 10)$y) %>% t()
 
  dat <- data.frame(pos = seq_len(ncol(cov_mat)),
                    height = colMeans(cov_mat, na.rm = T) )
  
  ggplot(dat, aes(x = pos, y = height)) +
    geom_vline(xintercept = ncol(cov_mat) %/% 2 + 1, color = "grey80", lwd = 1) +
    geom_line(lwd = 1, color = line.color) +
    scale_x_continuous(breaks = c(1, ncol(cov_mat) %/% 2 + 1, ncol(cov_mat)),
                       labels = c(x_left_lab, x_mid_lab, x_right_lab)) +
    xlab(xlab) + ylab(ylab) + ggtitle(.title) +
    theme_setting
}


# plot
g2.1 <- plot_coverage_profile(log1p(trim_quantile(cov_TTseq)[, 3:27] ), line.color = colors_9[8],
                              xlab = "Termination site [kb]", ylab = "log mean coverage",
                              .title = "TT-seq nascent mRNA (n=10447)")

g2.2 <- plot_coverage_profile(log1p((cov_Pol2S5p[[1]])[, 3:27] ), line.color = colors_9[2],
                      xlab = "Termination site [kb]", ylab = "log mean coverage",
                      .title = "Pol II-S5p")

g2.3 <- plot_coverage_profile((log1p(cov_TTseq) - log1p(cov_Pol2S5p[[1]] ))[, 3:27], line.color = colors_9[3],
                              xlab = "Termination site [kb]", ylab = "log TT-seq / Pol II-S5p",
                              .title = "Est. speed (a.u.)")

g2.4 <- plot_coverage_profile(log1p(trim_quantile(cov_FAIRE)[, 3:27]), line.color = colors_9[4],
                              xlab = "Termination site [kb]", ylab = "log mean coverage",
                              .title = "FAIRE")

# daRNA transcription collision 
daRNA.gr <- TU.DE.mm9[TU.DE.mm9$location %in% c("daRNA", "intergenic")]
cov_da <- intervalCoverage(list(promoters(daRNA.gr, upstream = 100, downstream = 100)),
                           intervals = terSite_SL + 1999, out_width = 25)

g2.5 <- plot_coverage_profile(cov_da[[1]][[2]][rowSums(cov_da[[1]][[2]]) > 0, ],
                              line.color = colors_9[5], x_left_lab = "-2", x_right_lab = "2",
                              xlab = "Termination site [kb]", ylab = "Occurence",
                              .title = "daRNA TSS (n=714)")


ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5, nrow = 1),
       filename = paste0("FigS5_terSite_coverages.png"), 
       path = "../figS5/figs/",
       device = "png", width = 20, height = 4)
