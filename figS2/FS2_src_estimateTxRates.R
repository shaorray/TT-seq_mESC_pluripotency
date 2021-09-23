#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
## predict copy number per cell
weightTransform <- function( W ) {
  out_W <- 2^W - 0.001
  out_W[out_W < 0] <- 0
  out_W
}

# Constants for RNA copy number transformation
# 1) correct spike-in weight due to the different mole number design to adapt to 
# the final mixed concentration 
#     before      now
# Sp2      1        1
# Sp4      1      0.1
# Sp5      1        1
# Sp8      1      0.1
# Sp9      1     0.01
# Sp12     1     0.01
# sum      6     2.22
norm_factor <- 6 / 2.22 

# 2) Transform 0.4 ng spike-in RNA mixture per million cell to molecular number per cell:
# ( 0.4 * 1e-9 (g) / 322 (nt molecular mass) * 6.02e+23 ) / 1000 (to RPK average length) / 1e+6 (million cells) / norm_factor 
transform_facter <- 0.4e-9 / 322 / 1000 * 6.02e+23 * 1e-6 / norm_factor

## predict labeled rate and half-life
PredictRate <- function(beta, LRNA, FRNA) {
  rate_term <- cbind(1, FRNA, LRNA) %*% beta
  rate <- 2^(rate_term) - 0.001
  rate[is.na(rate) | rate < 0] <- 0
  rate
}

fit_weight <- readRDS("data/fit_weight_lm.RData")
beta_r <- coef(readRDS("data/labeled_rate_lm_fit.RData"))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# geomatric mean by row
geoMeans <- function(x) exp(log(matrixStats::rowProds(x)) / ncol(x))

sampleMatMeans <- function(sample_counts, sample_names, FUN) {
  # cbind(sample_counts[, 1:2] %>% FUN,
  #       sample_counts[, 3],
  #       sample_counts[, 4:5] %>% FUN,
  #       sample_counts[, 6:7] %>% FUN,
  #       sample_counts[, 8:9] %>% FUN
  # )
  sapply(paste0(sample_names, "_rep"), 
         function(x) {
           tmp <- sample_counts[, grep(x, colnames(sample_counts))]
           if (!is.null(dim(tmp))) {
             tmp %>% FUN
           } else {
             tmp
           }
         }
    )
}

# ------------------------------------------------------------------------
get_label_rate_copy_mat <- function(L_mat, F_mat) {
  # melt data columns
  sample_counts <- NULL
  for(i in colnames(L_mat)){
    sample_counts <- rbind(sample_counts,
                           data.frame(FRNA = F_mat[, i],
                                      LRNA = L_mat[, i],
                                      Sample = i,
                                      Id = rownames(F_mat)))
  }
  
  ## predict copy number
  sample_counts$Copy <- predict(fit_weight,
                                data.frame(FRNA = sample_counts$FRNA,
                                           LRNA = sample_counts$LRNA)
  ) %>% weightTransform * transform_facter
  
  ## calculate labeled rate and half-life
  sample_counts$R <- with( sample_counts,
                           PredictRate( beta = beta_r,
                                        FRNA = log2(FRNA),
                                        LRNA = log2(LRNA) )
  )
  sample_counts$R <- with(sample_counts, ifelse(FRNA == 0 | LRNA == 0, NA, R)) # remove low quality estimation
  
  # reorder samples
  rate_mat <- matrix(sample_counts$R, ncol = ncol(L_mat))
  
  sample_names <- unique(gsub("\\_rep.*", "\\1", sample_counts$Sample))
  num_genes <- length(unique(sample_counts$Id))
  
  sample_counts_merged <- 
    data.frame(Sample = rep(sample_names, each = num_genes),
               gene_id = rep(unique(sample_counts$Id), length(sample_names)),
               
               FRNA = c(sampleMatMeans(F_mat, sample_names, rowMeans)),
               LRNA = c(sampleMatMeans(L_mat, sample_names, rowMeans)),
               
               Copy = c(sampleMatMeans(`colnames<-`(matrix(sample_counts$Copy, ncol = ncol(L_mat)), colnames(L_mat)), 
                                       sample_names, rowMeans) ),
               Copy_sd = c(sampleMatMeans(`colnames<-`(matrix(sample_counts$Copy, ncol = ncol(L_mat)), colnames(L_mat)), 
                                          sample_names, matrixStats::rowSds ) ),
               
               Labeled_rate = c(sampleMatMeans(`colnames<-`(rate_mat, colnames(L_mat)), sample_names, geoMeans) ),
               Labeled_rate_sd = c(sampleMatMeans(`colnames<-`(matrix(sample_counts$R, ncol = ncol(L_mat)), colnames(L_mat)),
                                                  sample_names, matrixStats::rowSds ) )
  )
  sample_counts_merged$Half_life <- c( log(2) / log(1 - sample_counts_merged$Labeled_rate) * (-5) )
  
  sample_counts_merged
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GENCODE transcripts half-life and copy
tx_L_mat <- readRDS("../data/tx_L_mat.RData") 
tx_F_mat <- readRDS("../data/tx_F_mat.RData")
sample_Tx_counts_merged <-
  get_label_rate_copy_mat(L_mat = tx_L_mat[, !grepl("rep3|7d", colnames(tx_L_mat))],
                          F_mat = tx_F_mat[, !grepl("rep3|7d", colnames(tx_F_mat))])

saveRDS(sample_Tx_counts_merged, 'data/sample_Tx_counts_Rates_combined2.RData')
write.table(sample_Tx_counts_merged, file = '../data/Sample_RNA_counts_rates_mat2.txt', quote = F, sep = '\t', row.names = F)
# write.table(unique(sample_Tx_counts_merged$gene_id[with(sample_Tx_counts_merged, FRNA * LRNA > 0)]),
#             file = '../data/active_es_genes.txt',
#             quote = F, sep = '\t', row.names = F, col.names = F)

cat(paste(table(sample_Tx_counts_merged$Sample)[1], "TXs median half-life by sample: \n")) # about 1 h
print(aggregate(sample_Tx_counts_merged$Half_life,
                by = list(sample_Tx_counts_merged$Sample),
                median, na.rm = T))
cat("\n")

# ---------------------------------------------------------------------------------------------
# ncRNA reads counting with mapped bam files gives similar results with kallisto pseudo reads
# counts on the GENECODE TX references and ncRNA TU annotated sequences

# Annotated non-coding transcripts half-life and copy
tu_L_mat <- readRDS("../fig1/data/tu_LRNA_RPK_norm.RData") # kallisto tpm, spike-in normalized
tu_F_mat <- readRDS("../fig1/data/tu_FRNA_RPK_norm.RData")
sample_nc_counts_merged <- 
  get_label_rate_copy_mat(L_mat = tu_L_mat[, !grepl("rep3|7d", colnames(tu_L_mat))],
                          F_mat = tu_F_mat[, !grepl("rep3|7d", colnames(tu_F_mat))])

saveRDS(sample_nc_counts_merged, 'data/sample_nc_counts_Rates_combined2.RData')

cat(paste(table(sample_nc_counts_merged$Sample)[1], "ncRNAs median half-life by sample: \n")) # about 30 min
print(aggregate(sample_nc_counts_merged$Half_life,
                by = list(sample_nc_counts_merged$Sample),
                median, na.rm = T))
cat("\n")

# -------------------------------------------------------------------------------------- #
# estimate cross-contamination rates
if (F) {
  # FRNA (X), LRNA (Y), labeled rate (R), cross-contamination rate (C)
  SampleSpCounts <- readRDS("../data/SampleSpikeCounts.RData")
  #1 scale fully labeled FRNA and LRNA (Sp2,4) in all samples
  lm1 <- lm(log2(LRNA) ~ I(log2(R*FRNA)), SampleSpCounts[SampleSpCounts$R > 0, ])
  
  #2 predict Sp9,12 theoretical labeled RNA if fully labeled
  LRNA1 <- predict(object = lm1,
                   newdata = data.frame(FRNA = SampleSpCounts[SampleSpCounts$R == 0, "FRNA"],
                                        R = 1))
  X_dat <- cbind(SampleSpCounts[SampleSpCounts$R == 0, c("Sample", "SpikeIns")],
                 X_contamination = SampleSpCounts[SampleSpCounts$R == 0, "LRNA"] / 2^LRNA1)
  
  ggplot(X_dat, aes(x = Sample, y = X_contamination * 100, 
                    group = Sample, color = SpikeIns)) +
    geom_line(color = "grey50") +
    geom_point(size = 4) +
    scale_color_manual(name = "UnlabeledSpikeins", values = colors_9) +
    xlab("") + ylab("Cross contamination %") + 
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))
  
  ggsave(filename = "FigS2_dot_sample_cross_contamination.png",
         device = "png", path = "../figS2/figs", width = 5, height = 4)
}
