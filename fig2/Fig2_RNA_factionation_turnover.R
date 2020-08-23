# RNA turnover validation with subcellular fractionation
# Rui Shao
# Feb 2020
# ----------------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
require(gridExtra)

source("F2_src_cell_fractionation_RNA.R") # rerun if data below not exist
# ----------------------------------------------------------------------------------------
fit_data <- readRDS("data/fractionation_fit_data.RData")
fit_data$half_life_k = log(2) / fit_data$k
fit_data$half_life_r = log(2) / fit_data$r
fit_data = fit_data[fit_data$half_life_k < 1000 & 
                      fit_data$half_life_r < 1000 &
                      fit_data$half_life_k > 30 &
                      fit_data$half_life_r > 5, ] # 4836

# measured total RNA level ralative to theoretical steady state
rho = (1 - exp(-5 * fit_data$k)) / (1 - exp(-5 * fit_data$r)) 
# if let nascent RNA decay equals to 0, results are very similar
# rho = 5 * fit_data$k / (1 - exp(-5 * fit_data$r)) 
names(rho) = rownames(fit_data)

# half synthsis replacement time: log(1 + rho) / fit_data$k

# show the correlation of two methods
ggplot(fit_data, aes(x = log10(half_life_r), y = log10(half_life_k))) +
  geom_hex(bins = 50) +
  scale_x_continuous(name="Synthesis Half-life (min)",
                     breaks=1:3, labels=10^(1:3)) +
  scale_y_continuous(name="Decay Half-life (min)",
                     breaks=c(1.5, 2, 2.5, 3), labels=round(10^c(1.5, 2, 2.5, 3), 0)) +
  geom_text(x = 2.7, y = 1.7, label = paste("r =", round(cor(log10(fit_data$half_life_k), log10(fit_data$half_life_r)),3) )) +
  geom_text(x = 2.7, y = 1.6, label = paste("n =", nrow(fit_data) )) +
  scale_fill_gradientn(colours = rev(colors_n)) +
  theme_setting +
  theme(legend.position = "none")

ggsave(filename = "Fig2_label_rate_decay_rate_half_life_comparison.pdf",
       path = "figs", device = 'pdf', width = 4.5, height = 4)

# steady state relative level
rho <- sort(rho)
high_GOs <- names(rho)[rho > 1.2]
low_GOs <- names(rho)[rho < 0.08]
data.frame(F_x = seq_along(rho) / length(rho),
           rho = log10(rho),
           colors = ifelse(rho > 1.2, "high", ifelse(rho < 0.08, "low", "median"))) %>%
  ggplot(aes(x = rho, y = F_x, color = colors) ) + 
  geom_line(lwd = 2) + 
  scale_color_manual(values = c("red", "blue", "darkgrey")) +
  scale_x_continuous(breaks = c(-2, -1, log10(median(rho)), 0, 1),
                     labels = c(10^c(-2, -1), round(median(rho), 2), 10^(0:1)),
                     limits = c(-2, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.91, 1), 
                     labels = c(0, 0.25, 0.5, 0.75, 0.91, 1)) +
  geom_segment(aes(x = -2, y = 0.91, xend = 0, yend = 0.91),
               col = "grey", lty = 2) +
  geom_segment(aes(x = -2, y = 0.5, xend = -0.44, yend = 0.5),
               col = "grey", lty = 2) +
  geom_segment(aes(x = -0.44, y = 0, xend = -0.44, yend = .5),
               col = "grey", lty = 2) +
  xlab("\nTotal RNA / Steady state level") + ylab("ECDF") +
  theme_setting + 
  theme(legend.position = "none", axis.title.x = element_text(size = 13))

ggsave(filename = "Fig2_total_RNA_steady_state_ratio.pdf",
       path = "figs", device = 'pdf', width = 3, height = 4)


mod_dat <- data.frame(x = seq(0,2,0.002), 
                      y = 1 - exp(-seq(0,2,0.002)),
                      colors = as.numeric(table(cut(rho, 1 - exp(-seq(0,2.002,0.002))))) )
ggplot(mod_dat, aes(x = x, y = y, color = colors)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
               col = "grey", lty = 1) +
  geom_point(pch = 2, cex = 0.8) +
  scale_color_gradientn(colours = (c("#d6e7ff", "#5297f7", "#000000"))) + 
  scale_x_continuous(name = "Time (t)",
                     breaks = c(0, log(2)),
                     labels = c(0, "log(2) / λ" ),
                     limits = c(0, 2)) +
  scale_y_continuous(name = "Total RNA (Y)",
                     breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0, 1.05)) +
  geom_segment(aes(x = log(2), y = 0, xend = log(2), yend = 0.5),
               col = "grey", lty = 2) +
  geom_segment(aes(x = 0, y = 0.5, xend = log(2), yend = 0.5),
               col = "grey", lty = 2) +
  geom_segment(aes(x = 0, y = 1, xend = 2, yend = 1),
               col = "grey", lty = 2) +
  geom_segment(aes(x = 0.36, y = 0, xend = 0.36, yend = 0.3023237),
               col = "grey", lty = 1, lwd = 0.4) +
  geom_segment(aes(x = 0.33, y = 0, xend = 0.33, yend = 0.2810763),
               col = "grey", lty = 1, lwd = 0.4) +
  geom_segment(aes(x = 0, y = 0.3023237, xend = 0.36, yend = 0.3023237),
               col = "grey", lty = 1, lwd = 0.4) +
  geom_segment(aes(x = 0, y = 0.2810763, xend = 0.33, yend = 0.2810763),
               col = "grey", lty = 1, lwd = 0.4) +
  geom_text(x = 0, y = 1.04, label = "Steady state", cex = 4.5, hjust = 0, col = "black") +
  geom_text(x = 1, y = 0.97, label = "Synthesis X = μt", cex = 4.5, hjust = 0, col = "black") +
  geom_text(x = 0.4, y = 0.1, label = "if Δt << log(2) / λ", cex = 4.5, hjust = 0, col = "black") +
  geom_text(x = 0, y = 0.35, label = "ΔY→0", cex = 4.5, hjust = 0, col = "blue2") +
  geom_text(x = 0.5, y = 0.35, label = "X(Δt)→μ*Δt", cex = 4.5, hjust = 0, col = "blue2") +
  theme_setting +
  theme(legend.position = "none")
  
ggsave(filename = "FigS2_schematic_synthesis_steady_state.png",
       path = "../figS2/figs", device = 'png', width = 5, height = 4)


# subcellular fraction class ~ labeled rate
tx_label_rate <- readRDS("data/SL_tx_label_rate.RData")
gene_attributes <- readRDS("data/Fractionation_gene_attributes.RData")
gene_attributes_label_rate <- gene_attributes[match(tx_label_rate$gene_id, gene_attributes$gene_id), ]

gene_attributes_label_rate$half_life <-
  tx_label_rate$Labeled_rate[match(gene_attributes_label_rate$gene_id, tx_label_rate$gene_id) ]
gene_attributes_label_rate$half_life <- log(2) / log(1 - gene_attributes_label_rate$half_life) * (-5) 

data.frame(Log_labeled_rate = log10(gene_attributes_label_rate$half_life),
           RNA_type = gene_attributes_label_rate$RNA_type) %>%
  ggplot(aes(x = Log_labeled_rate)) +
  geom_histogram(aes(fill = RNA_type), position = "dodge", bins = 50) +
  xlab('log(labeled rate)') +
  ylab('Genes') +
  scale_x_continuous(name="Half-life (min)",
                     breaks=0:3, labels=10^(0:3), limits = c(0, 3)) +
  geom_text(x = .5, y = 520, label = paste("n =", nrow(gene_attributes_label_rate) )) +
  scale_fill_manual( values = colors_20[c(2, 6)] ) +
  theme_setting ->
  plot1

data.frame(Log_labeled_rate = log10(gene_attributes_label_rate$half_life),
           Subcellular = gene_attributes_label_rate$subcellular_location) %>%
  ggplot(aes(x = Log_labeled_rate)) +
  geom_histogram(aes(fill = Subcellular), position = "dodge", bins = 50) +
  xlab('log(labeled rate)') +
  ylab('Genes') +
  scale_x_continuous(name="Half-life (min)",
                     breaks=0:3, labels=10^(0:3), limits = c(0, 3)) +
  geom_text(x = .5, y = 400, label = paste("n =", nrow(gene_attributes_label_rate) )) +
  scale_fill_manual( values = colors_20[c(12,20,1)] ) +
  theme_setting ->
  plot2

grid.arrange(plot1, plot2, ncol=1)
# ggsave(grid.arrange(plot1, plot2, ncol=1),
#        filename = "Fig2_label_rate_RNA_type_location_histogram.pdf",
#        path = "figs", device = 'pdf', width = 5, height = 4 )
ggsave(grid.arrange(plot1, plot2, ncol=1),
       filename = "Fig2_label_rate_RNA_type_location_histogram.png",
       path = "figs", device = 'png', width = 5, height = 4 )

# RNA labeled rate regression test ---------------------------------------------------------------
fit_label_rate <- lm( log(half_life_r) ~., data = fit_data[, c(1:6, 10)])
fit_decay_rate <- lm( log(half_life_k) ~., data = fit_data[, c(1:6, 9)])

summary(fit_label_rate)
summary(fit_decay_rate)

cbind("label_rate" = fit_label_rate$coefficients,
      "decay_rate" = fit_decay_rate$coefficients)

# plot coefficients ---------------------------------------------------------------
data.frame(Responses = c(rep("Synthesis Half-life", 5), rep("Decay Half-life", 5)),
           Coefficients = c(fit_label_rate$coefficients[-c(1, 7)], fit_decay_rate$coefficients[-c(1, 7)]),
           Subcellular_enrichment = rep(names(fit_decay_rate$coefficients[-c(1, 7)]), 2)) %>%
  ggplot(aes(x = Subcellular_enrichment, y = Coefficients, fill = Responses)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8, ) +
  scale_fill_manual( values = colors_20[c(1,2)] ) +
  xlab('') + ylab("Explanatory Variables Coefficients") +
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme_setting
ggsave(filename = "FigS2_label_rate_decay_rate_fractionation_regression.pdf",
       path = "../figS2/figs", device = 'pdf', width = 5, height = 4)
