# Descriptive analysis of TU annotations

# TTseq LRNA bam files processed with TU filter
# with setting:
#   Merge feature: protein_coding, licnRNA
#   STAN method: log poisson
#   Bin: 200 bp
#   Combine: exon

# Rui Shao
# 2020 Feb

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

trans <- function(x){pmin(x, 2) + 0.005*pmax(x-2,0)}
yticks <- c(0, 0.5, 1, 1.5, 2, 100)

# pdf("figs/TU reads abundance.pdf", 6, 5)
g1 <- ggplot(nascent_reads, aes(x = location, y = trans(Abundance), fill = Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("Transcription Abundance %") +
  scale_fill_manual(values = colors_20[c(13, 13, 13, 2, 2, 20, 20, 7, 7)]) +
  geom_rect(aes(xmin=0, xmax=6, ymin=2, ymax=2.2), fill="white") +
  scale_y_continuous(limits=c(0, NA), breaks=trans(yticks), labels=yticks) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
ggsave(g1, filename = "Fig1_TU_reads_abundance.pdf", path = "figs",
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
TU.list <- readRDS("data/TU.list.RData")
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

g3 <- ggplot(jaccard_list[["asRNA"]],
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
