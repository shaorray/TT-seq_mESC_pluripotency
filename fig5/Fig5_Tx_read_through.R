setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../util/utils.R')
source('../util/getCoverage.R')

# source("F5_src_last_exon.R")
# source("F5_src_termination_window_coverage.R")


# find termination sites using termination window (-5kb, 15kb) expanding the last exons
# include exon coverage 
# load data -----------------------------------------------------------------------------------
terWindow <- importRanges("data/terWindow.mm9.gff3")
gene.last.exon.gr <- importRanges("data/gene.last.exon.mm9.gff3")

TTseq.terWin.cov_norm <- readRDS("data/TTseq.terWin.cov_norm.RData")
Pol2S5p.MINUTE.terWin.cov <- readRDS("data/Pol2S5p.MINUTE.terWin.cov.RData")

read_thru_table <- readRDS("data/read_thru_table.RData")
# tx termination sites  -----------------------------------------------------------------------------------------------------
if (!file.exists("data/read_thru_table.RData")) {
  nascent.terWin.cov <- readRDS("data/nascent.terWin.cov_norm.RData")
  pol2.terWin.cov <- readRDS("data/pol2.terWin.cov_norm.RData")
  
  get_ter_site <- function(x) {
    if (any(is.na(scale(x)))) return(0)
    which.max(cumsum(scale(log1p(x))))
  }
  
  TTseq.ter.sites <- Reduce(cbind, lapply(TTseq.terWin.cov_norm, function(sample.cov) unlist(sapply(sample.cov, get_ter_site))))
  colnames(TTseq.ter.sites) <- names(TTseq.terWin.cov_norm)
  
  Pol2S5p.MINUTE.ter.sites <- Reduce(cbind,lapply(Pol2S5p.MINUTE.terWin.cov, function(sample.cov) unlist(sapply(sample.cov, get_ter_site))))
  colnames(Pol2S5p.MINUTE.ter.sites) <- names(Pol2S5p.MINUTE.terWin.cov)
  
  nascent.ter.sites <- Reduce(cbind, lapply(nascent.terWin.cov, function(sample.cov) unlist(sapply(sample.cov, get_ter_site))))
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

# plot a dropping example with Hsp90ab1  ---------------------------------------------------
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
  xlab("Termination window (Kb)") + ylab("Density difference") +
  theme(strip.text.y = element_text(size = 14),
        plot.margin = margin(t = 0, r = 0.8, b = 0, l = 0.15, unit = "cm"))

ggsave(plot = grid.arrange(g1, g2, nrow = 2, heights = c(9, 3)),
       filename = "Fig5_terWindow_cov_example.png", path = "figs",
       device = "png", width = 5, height = 7)

g3 <- data.frame(x = read_thru_table$LRNA_2i / 1e3,
           y = read_thru_table$LRNA_SL / 1e3
           ) %>%
  ggplot(aes(x = x, y = y)) +
  geom_hex(bins=45) +
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="2i read through (kb)") +
  scale_y_continuous(name="SL read through (kb)") +
  theme(axis.title=element_text(size=16, face="bold"),
        legend.position = "none")

g4 <- data.frame(x = read_thru_table$LRNA_mTORi / 1e3,
                 y = read_thru_table$LRNA_SL / 1e3
) %>%
  ggplot(aes(x = x, y = y)) +
  geom_hex(bins=45) +
  theme_setting +
  scale_fill_viridis(option = "A", direction = -1) +
  scale_x_continuous(name="mTORi read through (kb)") +
  scale_y_continuous(name="SL read through (kb)") +
  theme(axis.title=element_text(size=16, face="bold"),
        legend.position = "none")

ggsave(plot = grid.arrange(g3, g4, nrow = 1),
       filename = "Fig5_terWindow_SL_2i_mTORi.png", path = "figs",
       device = "png", width = 7, height = 3.5)

# elongation speed ~ read through length  ---------------------------------------------------
elongation_table <- read.table('../data/elongation_table.txt')
read_thru_table <- readRDS("data/read_thru_table.RData")

gene_ids <- intersect.Vector(read_thru_table$gene_id, rownames(elongation_table))
mtch <- match(gene_ids, read_thru_table$gene_id)
elongation_table <- elongation_table[gene_ids, ] %>%
  add_column(read_through_SL = read_thru_table$LRNA_SL[mtch],
             read_through_2i = read_thru_table$LRNA_2i[mtch],
             read_through_mTORi = read_thru_table$LRNA_mTORi[mtch]) %>% 
  as.data.frame() %>%
  dplyr::filter(speed > 0.1 & speed < 12  & read_through_SL < 15000 & read_through_2i < 15000 & read_through_mTORi < 15000)

elongation_table$Speed_cls <- cut(log(elongation_table$speed),
                                  breaks = quantile(log(elongation_table$speed), seq(0, 1, 0.2)),
                                  labels = c("Very slow", "Slow", "Medium", "Fast", "Very fast") )
saveRDS(elongation_table, "data/elongation_table.RData")

g5 <- ggplot(elongation_table, aes(x = speed, y = read_through_SL / 1000 )) +
  # geom_point(color = add.alpha("black", 0.1), size = 0.7) +
  geom_hex(bins = 45) +
  geom_density2d(color = "black", size = 0.4, bins = 10) +
  # geom_smooth(method = 'glm', size = 0.6, color = "black") +
  scale_fill_gradientn(colours = add.alpha(rev(colors_n), 0.6)) +
  scale_x_log10() +
  geom_text(x = 0.85, y = 12.5, label = paste("r =", round(with(elongation_table, cor((speed), (read_through_SL))), 3) )) +
  geom_text(x = 0.85, y = 11.5, label = paste("p <", formatC(cor.test(log10(elongation_table$speed), elongation_table$read_through_SL)$p.val, format = "e", digits = 0) )) +
  geom_text(x = 0.85, y = 10.5, label = paste("n =", nrow(elongation_table)) ) +
  ylab("Termination window (Kb)") +
  xlab("\nElongation speed (Kb/min)\n") +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), legend.position = "none") +
  annotation_logticks(base = 10, sides = "bottom", scaled = T)


# RNA synthesis ~ read through length  ---------------------------------------------------
sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$Sample == "SL", ]
sample_Tx_counts <- sample_Tx_counts[match(rownames(elongation_table), sample_Tx_counts$gene_id), c("Copy", "Labeled_rate")]

dat_synthesis_RT <- data.frame(read_through_SL = elongation_table$read_through_SL,
                               synthesis = sample_Tx_counts %>% log10() %>% rowSums() %>% "-"(log10(5)) )
dat_synthesis_RT <- dat_synthesis_RT[!is.na(dat_synthesis_RT$synthesis) & !is.na(dat_synthesis_RT$read_through_SL), ] 

g6 <- ggplot(dat_synthesis_RT, aes(x = synthesis, y = read_through_SL / 1000 )) +
  # geom_point(color = add.alpha("black", 0.1), size = 0.7) +
  geom_hex(bins = 50) +
  geom_density2d(color = "black", size = 0.4, bins = 10) +
  scale_fill_gradientn(colours = add.alpha(rev(colors_n), 0.4)) +
  geom_text(x = -0.1, y = 12.5, label = paste("r =", round(with(dat_synthesis_RT , cor(10^(synthesis), (read_through_SL))), 3) )) +
  geom_text(x = -0.1, y = 11.5, label = paste("p <", formatC(cor.test(10^(dat_synthesis_RT$synthesis), dat_synthesis_RT$read_through_SL)$p.val, format = "e", digits = 0) )) +
  geom_text(x = -0.1, y = 10.5, label = paste("n =", nrow(dat_synthesis_RT)) ) +
  scale_x_continuous(breaks = (-2:0), labels = 10^(-2:0), limits = c(-2.2, 0.2)) +
  ylim(c(0.1, 13)) +
  xlab("\nNascent RNA (copy/min)\n") +
  ylab("Termination window (Kb)") +
  theme_setting +
  theme(axis.ticks.x = element_blank(), legend.position = "none") +
  annotation_logticks(base = 10, sides = "bottom", scaled = T)

ggsave(plot = grid.arrange(g5, g6, nrow = 2),
       filename = "Fig5_elongation_speed_synthesis_termination_window.png", path = "figs",
       device = "png", width = 3.5, height = 7)

# elongation speed class ~ ttseq coverage ---------------------------------------------------
names(gene.last.exon.gr) <- gene.last.exon.gr$gene_id
genes.intercect <- intersect.Vector(terWindow[width(terWindow) == 15000]$gene_id, intersect.Vector(gene.last.exon.gr$gene_id, rownames(elongation_table)))
elongation_table <- elongation_table[genes.intercect, ]
pr_genes_terWin <- gene.last.exon.gr[genes.intercect]

ttseq.cov.exonlast.list = readBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                         pattern = "LRNA.*SL_(rep1|rep2).*bam$", full.names = T),
                                  intervals = pr_genes_terWin,
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
  ylab("LRNA coverage") +
  scale_y_continuous(breaks = c(0.15, 1.15), labels = c("0.0", "1.0"),
                     limits = c(0, 1.15), expand = c(0,0)) +
  theme_setting +
  theme(axis.text=element_text(size=10, face = "plain"),
        axis.title=element_text(size=14,face="bold"),
        strip.text.y = element_text(size = 12, face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 0.5, lineend = "butt"),
        axis.ticks.y.left = element_line(size = 0.5),
        panel.spacing = unit(1.5, "lines"))
ggsave(filename = "Fig5_Elongation_speed_read_thru_LRNA.png", device = "png",
       width = 3.5, height = 7, path = "figs")

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
plot(density(read_thru_table$LRNA_SL), axes = T, lwd = 15, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()

png("figs/Fig5_Speed_coverage_2i.png", width = 1500, height = 500)
image(t(log1p(speed_cov_2i[order(read_thru_table$LRNA_2i %/% 100, decreasing = F), ])),
      xaxt = "none", yaxt = "none", col = colors, breaks = breaks)
axis(side=1, at=c(0, 0.333, 0.667, 1), labels=c(0, 5, 10, 15))
box()
dev.off()
png("figs/Fig5_Speed_terSite_density_2i.png", width = 1500, height = 500)
plot(density(read_thru_table$LRNA_2i), axes = T, lwd = 15,  main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()

png("figs/Fig5_Speed_coverage_mTORi.png", width = 1500, height = 500)
image(t(log1p(speed_cov_mTORi[order(read_thru_table$LRNA_mTORi %/% 100, decreasing = F), ])),
      xaxt = "none", yaxt = "none", col = colors, breaks = breaks)
axis(side=1, at=c(0, 0.333, 0.667, 1), labels=c(0, 5, 10, 15))
box()
dev.off()
png("figs/Fig5_Speed_terSite_density_mTORi.png", width = 1500, height = 500)
plot(density(read_thru_table$LRNA_mTORi), axes = T, lwd = 15, main = '', xaxt = "n", yaxt = 'n', xlab = '', ylab = '')
dev.off()


png("figs/Fig5_Speed_coverage_scale.png", width = 100, height = 500)
image(x=1, y=seq(0, max(log1p(speed_cov_SL)), length.out=100),
      z=matrix(seq(0, max(log1p(speed_cov_SL)), length.out=100), 1, 100),
      col = colors,
      ylab='', xlab='', xaxt='n',  axes=T)
box()
dev.off()

# extract speed falls --------------------------------------------------------------------------------
get_slopes <- function(cov, site) coef(lm(cov[seq_len(site)] ~ seq_len(site) ))[2]
get_time <- function(cov, site) sum(1/cov[seq_len(site)])
slopes_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
                    function(x) get_slopes(cov = log1p(speed_cov_SL[x, ]), 
                                           site = read_thru_table$LRNA_SL[x] %/% 100 + 10))
names(slopes_SL) <- read_thru_table$gene_id

time_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
                  function(x) get_time(cov = speed_cov_SL[x, ], 
                                         site = read_thru_table$LRNA_SL[x] %/% 100 + 10))
names(time_SL) <- read_thru_table$gene_id

read_thru_speed_SL <- sapply(seq_len(nrow(speed_cov_SL)), 
                             function(x) 
                               mean(speed_cov_SL[x, seq_len(ceiling(read_thru_table$LRNA_SL[x] / 100))]))
names(read_thru_speed_SL) <- read_thru_table$gene_id
read_thru_speed_2i <- sapply(seq_len(nrow(speed_cov_2i)), 
                             function(x) 
                               mean(speed_cov_2i[x, seq_len(ceiling(read_thru_table$LRNA_2i[x] / 100))]))
read_thru_speed_mTORi <- sapply(seq_len(nrow(speed_cov_mTORi)), 
                             function(x) 
                               mean(speed_cov_mTORi[x, seq_len(ceiling(read_thru_table$LRNA_mTORi[x] / 100))]))

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

# read through explaination with estimated elongation speed, gene length, tx length, 2mers -----------------
estimated_speed <- readRDS("../fig4/data/samples_estimated_speed.RData")
read_thru_table <- readRDS("data/read_thru_table.RData")
rownames(read_thru_table) <- read_thru_table$gene_id
genes.intercect <- intersect.Vector(read_thru_table$gene_id, rownames(estimated_speed))

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

gene.list = getBM(attributes = c("ensembl_gene_id", "transcript_length", "percentage_gene_gc_content"), 
                  mart = useMart("mmusculus_gene_ensembl", biomart = "ensembl"))
gene.list <- gene.list[!duplicated(gene.list$ensembl_gene_id), ]
rownames(gene.list) <- gene.list$ensembl_gene_id

dat_test = cbind("Synthesis rate" = sample_Tx_counts[genes.intercect, "synthesis"], # log scale
                 "Gene length" = log(width(gene.gr[genes.intercect])), 
                 "Transcript length" = log(gene.list[genes.intercect, ]$transcript_length),
                 "G.B. GC%" = gene.list[genes.intercect, ]$percentage_gene_gc_content,
                 "R.T. GC%" = trim_quantile(read_thru_GC[genes.intercect]),
                 "Est. G.B. speed" = estimated_speed[genes.intercect, 1], # log scale
                 "Est. R.T. speed" = log(read_thru_speed_SL[genes.intercect]) %>% trim_quantile(0.975),
                 "R.T. speed slope" = slopes_SL[genes.intercect] %>% trim_quantile(0.975))
r2 = multi_variance_explained(dat_test,  log1p(read_thru_table[genes.intercect, 2])) 
names(r2) = colnames(dat_test)
dat_r2 <- data.frame(feature = colnames(dat_test), r2 = r2, color = c(1:8))
dat_r2$feature <- factor(dat_r2$feature, levels = rev(colnames(dat_test)))

ggplot(dat_r2, aes(x = feature, y = r2, fill = feature)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) +
  coord_flip() +
  ylab(expression(~R^2)) +
  labs(x = "Read through explained") +
  geom_text(y = 0.3, x = 8, label = paste("Σ R2 =", round(sum(r2), 3)), cex = 3.5) +
  geom_text(aes(label = round(r2, 3)), hjust = -0.05, cex = 3.5) +
  scale_fill_manual(values = rev(colors_20[c(19, 20, 20, 8,9, 7, 6, 5, 4, 3)])) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), labels = c(0, 0.1, 0.2, 0.3), limits = c(-0.01, 0.4)) + 
  theme_setting +
  theme(legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.position = "none")
ggsave(filename = "Fig5_Read_through_explained.png", path = "figs",
       device = "png", width = 3.5, height = 4)

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

# estimated elongation speed changes ~ ttseq coverage, and read through changes ---------------------------------------------------
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

sf <- SizeFactorCal(cbind(rowSums(ttseq.cov.exonlast.SL[, 51:100]),
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




# -----------------------------------------------------------------------------------------
# supplementory
# count reads on gene body intervals
gene.gr <- readRDS("../fig4/data/gene.gr.RData")
P_R_gene_idx <- match(rownames(elongation_table), gene.gr$gene_id) # pause release gene of the elongation speed table from GRO-seq

gene.body.gr <- GRanges()
for (i in seq_along(gene.gr)) {
  gene.tmp = gene.gr[i]
  .strand = as.character(strand(gene.tmp))
  if (.strand == "+") {
    start(gene.tmp) = start(gene.tmp) + 500
  } else {
    end(gene.tmp) = end(gene.tmp) - 500
  }
  gene.body.gr <- c(gene.body.gr, gene.tmp)
}

L_sf = readRDS("../data/LRNA.sizefactor.RData")
L_sf = L_sf[grep("SL_|^2i_2d|mTORi_2d", names(L_sf))]
TT_seq_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                           pattern = "LRNA.*(SL_|_2i_2d|mTORi_2d).*bam$", full.names = T)
TT_seq_counts <- .countBam(bam_files = TT_seq_files,
                           intervals = gene.body.gr,
                           stranded = T, 
                           paired.end = "midpoint")
TT_seq_counts <- sweep(TT_seq_counts, 2, L_sf, "/")
TT_seq_counts <- data.frame("LRNA_SL" = rowMeans(TT_seq_counts[, 5:7]),
                            "LRNA_2i" = rowMeans(TT_seq_counts[, 1:2]),
                            "LRNA_mTORi" = rowMeans(TT_seq_counts[, 3:4]))

Pol2S5p_sf <- readRDS("../fig4/data/pol2s5p_input_size_factor.RData")[c(3,1,4)]
Pol2S5p_files <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                            pattern = "P1.*5p.*ALL.*bam$", full.names = T)[c(3,1,4)]
Pol2S5p_counts <- .countBam(bam_files = Pol2S5p_files,
                            intervals = gene.body.gr,
                            stranded = F, 
                            paired.end = "midpoint")
Pol2S5p_counts <- sweep(Pol2S5p_counts, 2, Pol2S5p_sf, "/")
colnames(Pol2S5p_counts) <- c("Pol2S5p_SL", "Pol2S5p_2i", "Pol2S5p_mTORi")

other_methods_counts <- .countBam(bam_files = c(list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Jesse_PROSeq_ControlClone2_Rep1/bam",
                                                           pattern = "out.bam$", full.names = T),
                                                list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Jesse_PROSeq_ControlClone2_Rep2/bam",
                                                           pattern = "out.bam$", full.names = T),
                                                list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Marta_ESC_1_PRO-Seq_ctrl/bam",
                                                           pattern = "out.bam$", full.names = T),
                                                list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Marta_ESC_2_PRO-Seq_ctrl/bam",
                                                           pattern = "out.bam$", full.names = T),
                                                "Ferrai_Pol_II_S5p" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2017_Ferrai_ChIP_RNAPIIS5p_ESC/bam",
                                                                                 pattern = "Ferrai.*bam$", full.names = T),
                                                "Ferrai_Pol_II_S7p" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2017_Ferrai_ChIP_RNAPIIS7p_ESC/bam",
                                                                                 pattern = "Ferrai.*bam$", full.names = T),
                                                "Flynn_Pol_II_DMSO" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2016_Flynn_ChIPseq_mES_Pol_II_DMSO_JQ1/",
                                                                                 pattern = "DMSO.*bam$", full.names = T),
                                                "Zhang_Pol_II_S2p" = "/mnt/0E471D453D8EE463/GEO_pol2/2018_Zhang_E14_Pol2S2P_Rep1/bam/2018_Zhang_E14_Pol2S2P.mm9.bam",
                                                "Zhang_Pol_II_S5p" = "/mnt/0E471D453D8EE463/GEO_pol2/2018_Zhang_E14_Pol2S5P_Rep1/bam/2018_Zhang_E14_Pol2S5P.mm9.bam"
                                                ),
                            intervals = gene.body.gr,
                            stranded = T, 
                            paired.end = "ignore")

other_methods_counts <- data.frame("Jesse_PROSeq" = rowMeans(other_methods_counts[, 1:2]),
                                   "Marta_PROSeq" = rowMeans(other_methods_counts[, 3:4]),
                                   "Ferrai_Pol2S5p" = other_methods_counts[,5],
                                   "Ferrai_Pol2S7p" = other_methods_counts[,6],
                                   "Flynn_Pol2" = rowMeans(other_methods_counts[,7:8]),
                                   "Zhang_Pol2S2p" = other_methods_counts[,9],
                                   "Zhang_Pol2S5p" = other_methods_counts[,10] )
# convert to RPK
TT_seq_counts <- sweep(TT_seq_counts, 1, width(gene.body.gr)/1000, "/")
Pol2S5p_counts <- sweep(Pol2S5p_counts, 1, width(gene.body.gr)/1000, "/")
other_methods_counts <- sweep(other_methods_counts, 1, width(gene.body.gr)/1000, "/")

# SL speed vs read-through
dat = data.frame(v_hat = log(TT_seq_counts$LRNA_SL[P_R_gene_idx] / Pol2S5p_counts$Pol2S5p_SL[P_R_gene_idx]),
           read_through = elongation_table$read_through_SL / 1e3) 
dat = dat[!is.infinite(dat$v_hat), ]
g5 <- ggplot(dat, aes(x = v_hat, y = read_through)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "A", direction = -1, alpha = 0.8) +
  geom_text(x = 3.7, y = 13, label = paste("r =", round(with(dat, cor(v_hat, read_through)), 3) )) +
  ylab("Termination window (Kb)\n") +
  xlab("\nEstimated speed (a.u.)") +
  xlim(c(-2, 5)) +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none")

# 2i speed vs read-through
dat = data.frame(v_hat = log(TT_seq_counts$LRNA_2i[P_R_gene_idx] / Pol2S5p_counts$Pol2S5p_2i[P_R_gene_idx]),
                 read_through = elongation_table$read_through_2i / 1e3) 
dat = dat[!is.infinite(dat$v_hat), ]
g6 <- ggplot(dat, aes(x = v_hat, y = read_through)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "A", direction = -1, alpha = 0.8) +
  geom_text(x = 3.7, y = 13, label = paste("r =", round(with(dat, cor(v_hat, read_through)), 3) )) +
  ylab("Termination window (Kb)\n") +
  xlab("\nEstimated speed (a.u.)") +
  xlim(c(-2, 5)) +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank(), 
        legend.position = "none")

# mTORi speed vs read-through
dat = data.frame(v_hat = log(TT_seq_counts$LRNA_mTORi[P_R_gene_idx] / Pol2S5p_counts$Pol2S5p_mTORi[P_R_gene_idx]),
                 read_through = elongation_table$read_through_mTORi / 1e3) 
dat = dat[!is.infinite(dat$v_hat), ]
g7 <- ggplot(dat, aes(x = v_hat, y = read_through)) +
  geom_hex(bins = 50) +
  scale_fill_viridis(option = "A", direction = -1, alpha = 0.8) +
  geom_text(x = 3.7, y = 13, label = paste("r =", round(with(dat, cor(v_hat, read_through)), 3) )) +
  ylab("Termination window (Kb)\n") +
  xlab("\nEstimated speed (a.u.)") +
  xlim(c(-2, 5)) +
  ylim(c(0.1, 13)) +
  theme_setting +
  theme(axis.ticks.x = element_blank())

png(filename = "../figS5/figs/FigS5_Speed_termination_SL_2i_mTORi.png",
    height = 280, width = 900)
grid.arrange(g5, g6, g7, nrow = 1, widths = c(4, 4, 4.8))
dev.off()

