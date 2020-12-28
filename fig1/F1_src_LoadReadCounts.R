#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# read in kallisto tx counts on gencode.vM17.annotation and spike-in RNAs
# and save reads counts for normalization
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

# load kallisto counts, mm10
filenames <- sort(list.files('../data/kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                     as = 'matrix', what = "tpm")
colnames(count_table) <- sampleNewName

# split to GENCODE transcripts, annotated TUs, spike-in RNAs
txRPK <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T)

tuRPK <- count_table[!grepl("^chrS|^ENS", rownames(count_table)), ]
tuRPK <- tuRPK[rowSums(tuRPK) > 0, ]

spRPK <- count_table[grep("^chrS", rownames(count_table)), ]

saveRDS(txRPK, "../data/txRPK_SL_2i.RData") # raw tmp without normalization
saveRDS(tuRPK, "../data/tuRPK_SL_2i.RData") # mm10 TU annotation kallisto tpm
saveRDS(spRPK, "../data/spRPK_SL_2i.RData")

# matrix for spike-in normalization
sp_F_mat <- spRPK[, grepl("FRNA", colnames(spRPK))]
colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
FRNA.sizefactor <- SizeFactorCal(sp_F_mat)

sp_L_mat <- spRPK[, grepl("LRNA",colnames(spRPK))]
colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins

saveRDS(FRNA.sizefactor, "../data/FRNA.sizefactor.RData")
saveRDS(LRNA.sizefactor, "../data/LRNA.sizefactor.RData")

if (F) { # spike-in size factors for read count normalization
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table) <- sampleNewName
  
  # norm to RPK
  eff_lengths <- SummarizedExperiment::readKallisto(paste0(filenames[1], "/abundance.tsv"),
                                                    as = 'matrix', what = "eff_length")
  count_table <- count_table / c(eff_lengths) * 1e3
  
  spRPK <- count_table[grep("^chrS", rownames(count_table)), ]
  
  sp_F_mat <- spRPK[, grepl("FRNA", colnames(spRPK))]
  colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
  FRNA.sizefactor <- SizeFactorCal(sp_F_mat)
  
  sp_L_mat <- spRPK[, grepl("LRNA",colnames(spRPK))]
  colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
  LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins
  
  saveRDS(FRNA.sizefactor, "../data/FRNA.sizefactor.RC.RData")
  saveRDS(LRNA.sizefactor, "../data/LRNA.sizefactor.RC.RData")
}

## spike-ins table, for labelel rate estimation
spikein_lens <- c("chrS2" = 1.023, "chrS4" = 1.033, "chrS5" = 1.042,
                  "chrS8" = 1.124, "chrS9" = 1.061, "chrS12" = 1.124) # convert RPK to abundance by multiplying spikein lengths
SampleSpCounts <- NULL
for(i in seq_len(ncol(sp_F_mat))){
  SampleSpCounts <- rbind(SampleSpCounts,
                          data.frame(FRNA = sp_F_mat[, i] / FRNA.sizefactor[i] * spikein_lens,
                                     LRNA = sp_L_mat[, i] / LRNA.sizefactor[i] * spikein_lens,
                                     Sample = colnames(sp_L_mat)[i],
                                     SpikeIns = rownames(sp_F_mat),
                                     W = c(1, 0.1, 1, 0.1, 1, 0.1),
                                     R = c(1, 1, 0.1, 0.1, 0, 0) ))
}
saveRDS(SampleSpCounts, "../data/SampleSpikeCounts.RData")

# normalise tx tpm with spike-in RNA size factor, for parameter estimation
tx_F_mat <- txRPK[, grep("FRNA", colnames(txRPK))]
colnames(tx_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tx_F_mat))
tx_F_mat <- sweep(tx_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tx_L_mat <- txRPK[, grep("LRNA", colnames(txRPK))]
colnames(tx_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tx_L_mat))
tx_L_mat <- sweep(tx_L_mat, 2, LRNA.sizefactor, '/')

saveRDS(tx_L_mat, "../data/tx_L_mat.RData")
saveRDS(tx_F_mat, "../data/tx_F_mat.RData")

# normalise non-coding TU tpm with spike-in RNA size factor
tu_F_mat <- tuRPK[, grep("FRNA", colnames(tuRPK))]
colnames(tu_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tu_F_mat))
tu_F_mat <- sweep(tu_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tu_L_mat <- tuRPK[, grep("LRNA", colnames(tuRPK))]
colnames(tu_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tu_L_mat))
tu_L_mat <- sweep(tu_L_mat, 2, LRNA.sizefactor, '/')

saveRDS(tu_F_mat, "../fig1/data/tu_FRNA_RPK_norm.RData")
saveRDS(tu_L_mat, "../fig1/data/tu_LRNA_RPK_norm.RData")

# process mm9 tx counts ------------------------------------------------------------------------------------------
if (T) {
  # count with combined reference GENCODE vM20 and ncRNA annotation
  filenames <- sort(list.files('../../figure_script/data/kallisto_output', full.names = T)) 
  sampleNewName <- gsub(".*/", "\\2", filenames)
  
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table) <- sampleNewName
  
  eff_lengths <- SummarizedExperiment::readKallisto(paste0(filenames[1], "/abundance.tsv"),
                                                    as = 'matrix', what = "eff_length")
  count_table <- count_table / c(eff_lengths) * 1e3
  
  txRPK <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
    keepOneTx(rowname_gene_id = T)
  txRPK <- txRPK[rowSums(txRPK) > 0, ]
  
  tuRPK <- count_table[!grepl("^chrS|^ENS", rownames(count_table)), ]
  tuRPK <- tuRPK[rowSums(tuRPK) > 0, ]
  
  spRPK <- count_table[grep("^chrS", rownames(count_table)), ]
  
  # matrix for spike-in normalization
  sp_F_mat <- spRPK[, grepl("FRNA", colnames(spRPK))]
  colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
  FRNA.sizefactor <- SizeFactorCal(sp_F_mat)
  
  sp_L_mat <- spRPK[, grepl("LRNA",colnames(spRPK))]
  colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
  LRNA.sizefactor <- SizeFactorCal(sp_L_mat[1:4, ]) # use only labeled spike-ins
  
  # normalise with spike-ins size factors
  tx_F_mat <- txRPK[, grepl("FRNA", colnames(txRPK))]
  colnames(tx_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tx_F_mat))
  tx_F_mat <- sweep(tx_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors
  
  tx_L_mat <- txRPK[, grepl("LRNA", colnames(txRPK))]
  colnames(tx_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tx_L_mat))
  tx_L_mat <- sweep(tx_L_mat, 2, LRNA.sizefactor, '/')
  
  tu_F_mat <- tuRPK[, grepl("FRNA", colnames(tuRPK))]
  colnames(tu_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tu_F_mat))
  tu_F_mat <- sweep(tu_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors
  
  tu_L_mat <- tuRPK[, grepl("LRNA",colnames(tuRPK))]
  colnames(tu_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tu_L_mat))
  tu_L_mat <- sweep(tu_L_mat, 2, LRNA.sizefactor, '/')
  
  # saveRDS(list("tx_F_mat" = tx_F_mat, "tx_L_mat" = tx_L_mat, 
  #              "tu_F_mat" = tu_F_mat, "tu_L_mat" = tu_L_mat), 
  #         "../fig1/data/tx_tu_RPK_norm.mm9.RData")
}

# ------------------------------------------------------------------------------------------
meanSampleCounts <- function(tx_mat)
{ # this function averaging replicates
  tmp_mat = NULL
  for( i in unique(colnames(tx_mat)) )
  {
    idx = colnames(tx_mat) == i

    if( sum(idx) > 1)
    {
      tmp_mat <- cbind(tmp_mat, rowMeans(tx_mat[, idx]))
    } else {
      tmp_mat <- cbind(tmp_mat, tx_mat[, idx])
    }
  }
  colnames(tmp_mat) <- unique(colnames(tx_mat))
  return(tmp_mat)
}
