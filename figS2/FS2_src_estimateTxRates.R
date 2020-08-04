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

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# geomatric mean by row
geoMeans <- function(x) exp( log(matrixStats::rowProds(x)) / ncol(x))

sampleMatMeans <- function( sample_counts, FUN ) {
  mat <- matrix(sample_counts, ncol = 9)
  cbind(sample_counts[, 1:2] %>% FUN,
        sample_counts[, 3],
        sample_counts[, 4:5] %>% FUN,
        sample_counts[, 6:7] %>% FUN,
        sample_counts[, 8:9] %>% FUN
  )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GENCODE transcripts half-life and copy
sample_Tx_counts <- NULL
for(i in colnames(tx_L_mat)[!grepl("rep3|SL2i", colnames(tx_L_mat))]){
  sample_Tx_counts <- rbind(sample_Tx_counts,
                            data.frame(FRNA = tx_F_mat[, i],
                                       LRNA = tx_L_mat[, i],
                                       Sample = i,
                                       Id = rownames(tx_F_mat)))
}

# ------------------------------------------------------------------------
# ncRNA reads counting with mapped bam files gives similar results with kallisto pseudo counts on
# the combined GeneCode TXs + ncRNA TU annotated sequences
## The following are counts from aligned reads

# # Append ncRNA read counts from bam featureCounts
# TU.TPM.mm9 <- readRDS("../fig1/data/TU.TPM.mm9.RData")
# colnames(mcols(TU.TPM.mm9)) <- gsub("\\.", "_",
#                                     gsub(".Aligned.sortedByCoord.out.bam", "\\2",
#                                          colnames(mcols(TU.TPM.mm9)))
#                                    )
# gene.idx <- which(TU.TPM.mm9$location == "protein_coding" & !is.na(TU.TPM.mm9$gene_id))
#
# # Normalise ncRNA according to kallisto tx TPM
# TU_L_mat = mcols(TU.TPM.mm9)[, grep("LRNA", colnames(mcols(TU.TPM.mm9)))]
# colnames(TU_L_mat) = gsub(".RNA_", "\\2", colnames(TU_L_mat))
# TU_L_mat = TU_L_mat[, colnames(tx_L_mat)]
# sf = cbind(tx_L_mat[match(TU.TPM.mm9$gene_id[gene.idx], rownames(tx_L_mat)), ],
#            TU_L_mat[gene.idx, ]) %>%
#   as.matrix() %>%
#   SizeFactorCal()
# TU_L_mat = t( t(as.matrix(TU_L_mat)) /
#                 sf[seq_len(length(sf) / 2) + length(sf) / 2] * sf[seq_len(length(sf) / 2)]
#             )
#
# TU_F_mat = mcols(TU.TPM.mm9)[, grep("FRNA", colnames(mcols(TU.TPM.mm9)))]
# colnames(TU_F_mat) = gsub(".RNA_", "\\2", colnames(TU_F_mat))
# TU_F_mat = TU_F_mat[, colnames(tx_L_mat)]
# sf = cbind(tx_F_mat[match(TU.TPM.mm9$gene_id[gene.idx], rownames(tx_F_mat)), ],
#            TU_F_mat[gene.idx, ]) %>%
#   as.matrix() %>%
#   SizeFactorCal()
# TU_F_mat = t( t(as.matrix(TU_F_mat)) /
#                 sf[seq_len(length(sf) / 2) + length(sf) / 2] * sf[seq_len(length(sf) / 2)]
# )
#
# ncRNA_counts <- NULL
# ncRNA.idx <- TU.TPM.mm9$location != "protein_coding"
# ncRNA.mm9 <- TU.TPM.mm9[ncRNA.idx]
# for(i in seq_len(ncol(TU_F_mat))){
#   ncRNA_counts <- rbind(ncRNA_counts, data.frame(FRNA = TU_F_mat[ncRNA.idx, i],
#                                                      LRNA = TU_L_mat[ncRNA.idx, i],
#                                                      Sample = colnames(TU_F_mat)[i],
#                                                      Id = paste(ncRNA.mm9$location,
#                                                                 seqnames(ncRNA.mm9),
#                                                                 start(ncRNA.mm9),
#                                                                 end(ncRNA.mm9),
#                                                                 strand(ncRNA.mm9),
#                                                                 sep = ";") ))
# }

# Annotated non-coding transcripts half-life and copy
tu_L_mat <- readRDS("../fig1/data/tu_LRNA_RPK_norm.rds") # kallisto tpm, spike-in normalized
tu_F_mat <- readRDS("../fig1/data/tu_FRNA_RPK_norm.rds")
ncRNA_counts <- NULL
for(i in colnames(tu_L_mat)[!grepl("rep3|SL2i", colnames(tu_L_mat))]){
  ncRNA_counts <- rbind(ncRNA_counts,
                        data.frame(FRNA = tu_F_mat[, i],
                                   LRNA = tu_L_mat[, i],
                                   Sample = i,
                                   Id = rownames(tu_L_mat) )) # on mm9 coordinate
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## predict copy number
fit_weight <- readRDS("data/fit_weight_lm.RData")

sample_Tx_counts$Copy <- predict(fit_weight,
                                 data.frame(FRNA = sample_Tx_counts$FRNA,
                                            LRNA = sample_Tx_counts$LRNA)
) %>% weightTransform * transform_facter

ncRNA_counts$Copy <- predict(fit_weight,
                             data.frame(FRNA = ncRNA_counts$FRNA,
                                        LRNA = ncRNA_counts$LRNA)
) %>% weightTransform * transform_facter

## calculate labeled rate and half-life
beta_r <- coef(readRDS("data/labeled_rate_lm_fit.RData"))
sample_Tx_counts$R <- with( sample_Tx_counts,
                            PredictRate( beta = beta_r,
                                         FRNA = log2(FRNA),
                                         LRNA = log2(LRNA) )
                          )
sample_Tx_counts$R <- with(sample_Tx_counts, ifelse(FRNA == 0 | LRNA == 0, NA, R)) # remove low quality estimation
sample_Tx_counts$Half_life <- log(2) / log(1- sample_Tx_counts$R) * (-5) # half-life minute

ncRNA_counts$R <- with( ncRNA_counts,
                        PredictRate( beta = beta_r,
                                     FRNA = log2(FRNA),
                                     LRNA = log2(LRNA) )
)
ncRNA_counts$R <- with(ncRNA_counts, ifelse(FRNA == 0 | LRNA == 0, NA, R))
ncRNA_counts$Half_life <- log(2) / log(1- ncRNA_counts$R) * (-5)

saveRDS(sample_Tx_counts, 'data/sample_Tx_counts_Rates.RData')
saveRDS(ncRNA_counts, 'data/sample_ncRNA_counts_Rates.RData')

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# reorder samples
rate_mat <- matrix(sample_Tx_counts$R, ncol = 9)

sample_names <- unique(gsub("\\_rep.*", "\\1", sample_Tx_counts$Sample))
num_genes <- length(unique(sample_Tx_counts$Id))

sample_Tx_counts_merged <- data.frame( Sample = rep(sample_names, each = num_genes),
                                       gene_id = rep(unique(sample_Tx_counts$Id), length(sample_names)),

                                       FRNA = c( sampleMatMeans(matrix(sample_Tx_counts$FRNA, ncol = 9), rowMeans) ),
                                       LRNA = c( sampleMatMeans(matrix(sample_Tx_counts$LRNA, ncol = 9), rowMeans) ),

                                       Copy = c( sampleMatMeans(matrix(sample_Tx_counts$Copy, ncol = 9), rowMeans) ),
                                       Copy_sd = c( sampleMatMeans( matrix(sample_Tx_counts$Copy, ncol = 9),
                                                                     matrixStats::rowSds ) ),

                                       Labeled_rate = c( sampleMatMeans(rate_mat, geoMeans) ),
                                       Labeled_rate_sd = c( sampleMatMeans( matrix(sample_Tx_counts$R, ncol = 9),
                                                                             matrixStats::rowSds ) )
                                      )

sample_Tx_counts_merged$Half_life <- c( log(2) / log(1 - sample_Tx_counts_merged$Labeled_rate) * (-5) )

saveRDS(sample_Tx_counts_merged, 'data/sample_Tx_counts_Rates_combined.RData')
write.table(sample_Tx_counts_merged, file = '../data/Sample_RNA_counts_rates_mat.txt', quote = F, sep = '\t', row.names = F)
write.table(unique(sample_Tx_counts_merged$gene_id[with(sample_Tx_counts_merged, FRNA * LRNA > 0)]),
            file = '../data/active_es_genes.txt',
            quote = F, sep = '\t', row.names = F, col.names = F)

cat(paste(table(sample_Tx_counts_merged$Sample)[1], "TXs median half-life by sample: \n")) # about 1 h
print(aggregate(sample_Tx_counts_merged$Half_life,
          by = list(sample_Tx_counts_merged$Sample),
          median, na.rm = T))
cat("\n")

# ncRNAs
nc_rate_mat <- matrix(ncRNA_counts$R, ncol = 9)

sample_names <- unique(gsub("\\_rep.*", "\\1", ncRNA_counts$Sample))
num_genes <- length(unique(ncRNA_counts$Id))

sample_nc_counts_merged <- data.frame( Sample = rep(sample_names, each = num_genes),
                                       gene_id = rep(unique(ncRNA_counts$Id), length(sample_names)),

                                       FRNA = c( sampleMatMeans(matrix(ncRNA_counts$FRNA, ncol = 9), rowMeans) ),
                                       LRNA = c( sampleMatMeans(matrix(ncRNA_counts$LRNA, ncol = 9), rowMeans) ),

                                       Copy = c( sampleMatMeans(matrix(ncRNA_counts$Copy, ncol = 9), rowMeans) ),
                                       Copy_Var = c( sampleMatMeans( matrix(ncRNA_counts$Copy, ncol = 9),
                                                                     matrixStats::rowVars ) ),

                                       Labeled_rate = c( sampleMatMeans(nc_rate_mat, geoMeans) ),
                                       Labeled_rate_Var = c( sampleMatMeans( matrix(ncRNA_counts$R, ncol = 9),
                                                                             matrixStats::rowVars ) )
)

sample_nc_counts_merged$Half_life <- c( log(2) / log(1 - sample_nc_counts_merged$Labeled_rate) * (-5) )

saveRDS(sample_nc_counts_merged, 'data/sample_nc_counts_Rates_combined.RData')

cat(paste(table(sample_nc_counts_merged$Sample)[1], "ncRNAs median half-life by sample: \n")) # about 30 min
print(aggregate(sample_nc_counts_merged$Half_life,
          by = list(sample_nc_counts_merged$Sample),
          median, na.rm = T))
cat("\n")
