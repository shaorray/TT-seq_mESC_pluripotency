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

cellCounts$Color <- factor(gsub("(.*)\\_rep.*", "\\1", cellCounts$Samples), c("SL", "2i_2d", "mTORi_1d", "mTORi_2d"))

ggplot(cellCounts, aes(x = CellNum / 145, y = L_sf - F_sf, color = Color)) + 
  geom_line(cex = 1.4, lty = 2) +
  geom_point(cex = 3) + 
  xlab("\nCell density (million / cm2)") +
  ylab("Turnover\n") +
  scale_colour_manual(values = colors_20[c(13, 2, 20, 7)]) +
  labs(color='Samples') +
  theme_setting
ggsave(filename = "Fig1.Turn_over_cell_number.pdf", 
       path = "../fig1/figs",
       device = "pdf", width = 5, height = 4 )

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
  geom_point() +
  geom_text(cex = 5, fontface = "bold", hjust = 0, nudge_x = 0.05) + 
  xlab("\nTotal RNA") +
  ylab("Nascent RNA\n") +
  xlim(c(-1, 1)) + 
  geom_hline(yintercept=0, cex = 0.2) + 
  geom_vline(xintercept=0, cex = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype="dotted") +
  geom_abline(intercept = 0, slope = -1, linetype="dotted") +
  scale_colour_manual(values = colors_20[c(2, 1, 20, 7, 13)]) +
  theme_setting +
  theme(legend.position = "none")
ggsave(filename = "Fig1.size_factors_distribution.pdf", path = "figs",
       device = "pdf", width = 4.5, height = 4.5 )
