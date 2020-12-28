# spike-in estimation RNA half-life and copy in absolute scales
# Rui Shao
# Nov 2019

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("EstimateTxRates.R")
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

# load data generated from "F1_src_LoadReadCounts.R"
tx_L_mat <- readRDS("../data/tx_L_mat.RData") 
tx_F_mat <- readRDS("../data/tx_F_mat.RData")
SampleSpCounts <- readRDS("../data/SampleSpikeCounts.RData")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate weight
fit_weight <- lm( log2(W + 0.001) ~ log2(FRNA) + log2(LRNA + 1), 
                  data = SampleSpCounts) # set lower limit of FRNA detection as 1/100 of sp2 copy
summary(fit_weight) # Adjusted R-squared:  0.9914
saveRDS(fit_weight, "data/fit_weight_lm.RData")

## plot and evaluate fitted weights
if (T) {
  # par(mfrow = c(2, 2))
  # plot(fit_weight)
  par(mfrow = c(1, 1))
  mat <- as.data.frame( matrix(0.001, nrow = 6, ncol = 3) ) # add near zero spike-ins pseudo number as the base line
  colnames(mat) = c("y", "FRNA", "LRNA")
  mat <- rbind( mat,
                data.frame(y = SampleSpCounts$W,
                           FRNA = SampleSpCounts$FRNA,
                           LRNA = SampleSpCounts$LRNA)
  )
  W <- predict(fit_weight, mat[, 2:3]) %>% weightTransform
  plot(x=mat[, 1], y= W, xlim = c(0, 1), ylim = c(0, 1.3),
       ylab = 'Predicted weights', xlab = 'Actual weights')
  abline(h = c(0, 0.1 , 1), col = 'red', lty = 2)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate labeled rate
# This formula follows the equation:
# Labeled RNA = Total RNA * labeled rate

fit_rate <- lm(log2(R + 0.001) ~  I(log2(FRNA + 1)) + I(log2(LRNA)), 
               data = SampleSpCounts[!SampleSpCounts$SpikeIns %in% c("chrS9", "chrS12"), ])
summary(fit_rate) # Adjusted R-squared: 0.9917

beta_r <- coef(fit_rate)
saveRDS(fit_rate, "data/labeled_rate_lm_fit.RData")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run the estimation of half-life and copy number
source("FS2_src_estimateTxRates.R")
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Test and plot fitted grid - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# leave-one-out, 10-fold cross validation of labeled rate
mat <- cbind(1,
             log2(SampleSpCounts$FRNA),
             log2(SampleSpCounts$LRNA),
             log2(SampleSpCounts$R + 0.001) )

cv_res <- NULL
group_residual <- NULL

for (i in 1:1000) {
  test_group <- sample(1:5, nrow(mat) %/% 5, T) + (seq_len(nrow(mat) %/% 5) - 1) * 5
  mat_train <- mat[-test_group, -4]
  beta_train <- solve(t(mat_train) %*% mat_train, t(mat_train) %*% mat[-test_group, 4])
  rate_pred <- PredictRate(beta_train, LRNA = mat[test_group, 3], FRNA = mat[test_group, 2])
  # Cross-validation results
  cv_res <- rbind( cv_res,
                   cbind(Labeled_rate = 2^(mat[test_group, 4]) - 0.001,
                         rate_pred = rate_pred ) )
  # Residual by sample
  group_residual <- rbind( group_residual, c(rate_pred - 2^(mat[test_group, 4]) + 0.001 ) )
}

par(mfrow = c(1, 1))

qqplot(cv_res[, 1], cv_res[, 2], xlab="Labeled rate", ylab="Labeled rate prediction", cex = 0.5)
abline(h=c(0, 0.1, 1), col='grey', lty=3)
colnames(group_residual) <- unique(SampleSpCounts$Sample)

par(mar=c(8,5,3,3))
boxplot(group_residual, outline=FALSE, ylim=c(-1, 1), ylab="Residuals", main="Predicted labeled rate by sample", las=2)
abline(h=0, col='grey', lty=3)

# saveRDS(group_residual, "data/rate_group_residual.RData")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## plot spike-ins half life grid, in supplement Fig S2

Fp <- seq(0, 17, by = 0.02) # FRNA part
Lp <- seq(0, 17, by = 0.02) # LRAN part

in_mat <- expand.grid(Fp, Lp)

out_mat <- matrix( PredictRate(beta = beta_r, FRNA = in_mat[, 2], LRNA = in_mat[, 1] ),
                   nrow = length(Fp) ) %>% t()

# out_mat <- apply(out_mat, 1, rev)
# out_mat <- apply(out_mat, 2, rev) 
out_mat[out_mat == 0 | out_mat >= 1] <- NA

# Convert to half-life and log10 transform for plotting
out_mat <- (log(2) / out_mat * 5) %>% log10() %>% log10()

if (T) {
  pdf("figs/FigS2_Spike_in_half_life_tile.pdf", 6, 6)
  par(mar=c(4,4,1,1))
  layout(mat=matrix(c(1,2), 1, 2),
         widths=c(2,0.6),
         heights=c(2,2), TRUE)

  cols = alpha(RColorBrewer::brewer.pal(9,'YlGnBu'), 0.9)

  image( out_mat, col=cols, axes=FALSE,
         xlab='log2(FRNA)', ylab='log2(LRNA)',
         main='Predicted grid of half-life')
  axis(1, at=c(0, 5, 10, 15) / 17, labels=c(0, 5, 10, 15), cex.axis=1)
  axis(2, at=c(0, 5, 10, 15) / 17, labels=c(0, 5, 10, 15), cex.axis=1)
  box()

  if (T) {
    par(new=TRUE)
    plot(0,xlim=c(0, 17), ylim=c(0, 17), type='n',axes=FALSE,ann=FALSE)
    points(log2(SampleSpCounts[, 1] + 1), log2(SampleSpCounts[, 2]),
           col=SampleSpCounts$SpikeIns)
    with(SampleSpCounts[1:6, ],
         text(x = log2(FRNA) + 0.5,y = log2(LRNA) + 0.8, adj = 0.1,
              labels = gsub("chr(.*)","\\1", SpikeIns), cex=1))

    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS2","chrS4"),],
                lm(log2(LRNA) ~ log2(FRNA)) ), lty= 2, lwd = 2)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS5","chrS8"),],
                lm(log2(LRNA) ~ log2(FRNA)) ), lty = 2, lwd = 2)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS2","chrS5", "chrS9"),],
                lm(log2(LRNA) ~ log2(FRNA))), lty = 2, lwd = 1.5)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS4","chrS8"),],
                lm(log2(LRNA) ~ log2(FRNA))), lty = 2, lwd = 1.5)
  }

  par(mar=c(4,1,1,4))
  image(x=1, y=seq(0, 10, length.out=200),
        z=matrix(seq(0, 10, length.out=200), 1, 200),
        col = cols,
        ylab='', xlab='Minute', xaxt='n',  axes=F)
  box()
  axis(4, at=seq(0, 10, 2.2222),
       labels = formatC( round(10^10^quantile(out_mat[!is.na(out_mat)],
                                              seq(0, 1, 0.2222)), digits = 2),
                         format = "e", 0),
       cex.axis=1, las=2)
  dev.off()
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## plot spike-ins copy number grid (per cell)

out_mat <- matrix( predict(fit_weight, data.frame(FRNA = 2^in_mat[, 2],
                                                  LRNA = 2^in_mat[, 1]) ) %>%
                     weightTransform,
                   ncol=length(Lp)) %>% t

## transform weight to molecular number
out_mat <- out_mat * transform_facter

if (T) {
  pdf("figs/FigS2_Spike_in_copy_number_tile.pdf", 6, 6)
  par(mar = c(4,4,1,1))
  layout(mat=matrix(c(1, 2), 1, 2),
         widths=c(2, 0.6),
         heights=c(2, 2), TRUE)

  cols = alpha(RColorBrewer::brewer.pal(9, 'YlOrBr'), 0.8)

  image( out_mat^0.2, col=cols, axes=FALSE, 
         xlab='log2(FRNA)', ylab='log2(LRNA)', 
         main='Predicted grid of weight')
  axis(1, at=c(0, 5, 10, 15) / 17, labels=c(0, 5, 10, 15), cex.axis=1)
  axis(2, at=c(0, 5, 10, 15) / 17, labels=c(0, 5, 10, 15), cex.axis=1)
  box()

  if (T) { # add spike_in counts
    par(new=TRUE)
    plot(0, xlim = c(0, 17), ylim=c(-1, 16), type='n', axes=FALSE, ann=FALSE)
    points(x = log2(SampleSpCounts[,1] + 1), y = log2(SampleSpCounts[,2] + 1),
           col = SampleSpCounts$SpikeIns)
    with(SampleSpCounts[1:6,],
         text(x = log2(FRNA) + 1, y = log2(LRNA) + 1, adj = 0.1,
              labels = gsub("chr(.*)","\\1", SpikeIns), cex=1))

    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS2","chrS4"),],
                lm(log2(LRNA) ~ log2(FRNA)) ), lty = 2, lwd = 2)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS5","chrS8"),],
                lm(log2(LRNA) ~ log2(FRNA)) ), lty = 2, lwd = 2)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS2","chrS5", "chrS9"),],
                lm(log2(LRNA) ~ log2(FRNA))), lty = 2, lwd = 1.5)
    abline(with(SampleSpCounts[SampleSpCounts$SpikeIns%in%c("chrS4","chrS8"),],
                lm(log2(LRNA) ~ log2(FRNA))), lty = 2, lwd = 1.5)
  }

  par(mar = c(4,1,1,4))

  image(x=1, y=seq(0, 10, length.out=300),
        z=matrix(seq(0, 10, length.out=300), 1, 300),
        col = cols,
        ylab='', xlab='Copy', xaxt='n',  axes=F)
  box()
  axis(4, at=seq(0, 10, 2.2222),
       labels = formatC( round( quantile( out_mat[!is.na(out_mat)]^0.2,
                                          seq(0, 1, 0.2222) ) %>%
                                  weightTransform, digits = 0),
                         format = "e", 0),
       cex.axis=1, las=2)
  dev.off()
}

# show spike-in half-life estimation versus conventional ratio half-life
sample_Tx_counts <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
sample_Tx_counts <- sample_Tx_counts[sample_Tx_counts$Sample == "SL", ]

g1 <- ggplot(sample_Tx_counts[with(sample_Tx_counts, Copy > 0.0001 & !is.na(Copy) & is.finite(Copy)), ], 
       aes(x = log10(FRNA), y = log10(Copy),
           color = get_dens(log10(FRNA), log10(Copy)))) +
  geom_point() +
  # scale_color_viridis_c(direction = -1) +
  scale_x_continuous(name = "log10 Total RNA RPK",
                     breaks = (-1:4),
                     labels = (-1:4),
                     limits = c(-1.5, 4.1)) +
  scale_y_continuous(name = "log10 Est. Copy",
                     breaks = (-3:5), 
                     labels = (-3:5), 
                     limits = c(-3, 4)) +
  annotate(geom = "text", 
           x = -1.5, y = 4,
           hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(with(sample_Tx_counts[with(sample_Tx_counts, Copy > 0.0001 & !is.na(Copy) & is.finite(Copy)), ], 
                                                   single_variance_explained(log10(FRNA), log10(Copy), is.cor = T)), 3)),
                          "\nn = ", sum(with(sample_Tx_counts, Copy > 0.0001 & !is.na(Copy) & is.finite(Copy))))) +
  theme_setting +
  theme(legend.position = "none")
  

g2 <- ggplot(sample_Tx_counts[with(sample_Tx_counts, !is.na(Half_life) & is.finite(Half_life)), ], 
       aes(x = log10(FRNA / LRNA), y = log10(Half_life),
           color = get_dens(log10(FRNA / LRNA), log10(Half_life)))) +
  geom_point() +
  # scale_color_viridis_c(direction = -1) +
  scale_x_continuous(name = "log10 Total RNA / Nascent RNA",
                     breaks = (-1:4),
                     labels = (-1:4),
                     limits = c(-2, 3)) +
  scale_y_continuous(name = "log10 Est. Half-life",
                     breaks = (-2:5),
                     labels = (-2:5),
                     limits = c(-1, 4)) +
  annotate(geom = "text", 
           x = -2, y = 4,
           hjust = "left", vjust = "top",
           label = paste0(paste("r = ", round(with(sample_Tx_counts[sample_Tx_counts$Copy > 0.0001, ], 
                                                   single_variance_explained(log10(FRNA / LRNA), log10(Half_life), is.cor = T)), 3)),
                          # "\np < ", formatC(with(elongation_table, 
                          #                        cor.test(log10(speed),
                          #                                 read_through_SL)$p.val), 
                          #                   format = "e", digits = 0),
                          "\nn = ", sum(with(sample_Tx_counts, !is.na(Half_life) & is.finite(Half_life))))) +
  theme_setting +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g1, g2, nrow = 1),
       filename = "FigS2_copy_half_life_estimation_versus_rpk_ratio.png", path = "figs",
       device = "png", width = 7, height = 3.5)
