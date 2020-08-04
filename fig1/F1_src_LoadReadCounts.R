# read in kallisto tx counts on gencode.vM17.annotation and spike-in RNAs
# and save reads counts for normalization
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

filenames <- sort(list.files('../data/kallisto_output/', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                     as = 'matrix', what = "tpm")
colnames(count_table) <- sampleNewName

txRPK <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T)
txRPK <- txRPK[apply(txRPK, 1, function(x) any(x != 0)), ]

tuRPK <- count_table[!grepl("^chrS|^ENS", rownames(count_table)), ]
tuRPK <- tuRPK[apply(tuRPK, 1, function(x) any(x != 0)), ]

spRPK <- count_table[grep("^chrS", rownames(count_table))[-6], ]

saveRDS(txRPK, "../data/txRPK_SL_2i.RData") # raw tmp without normalization
saveRDS(tuRPK, "../data/tuRPK_SL_2i.RData") # mm9 TU annotation tpm from kallisto
saveRDS(spRPK, "../data/spRPK_SL_2i.RData")

# matrix for spike-in normalization
sp_F_mat <- spRPK[,grepl("FRNA", colnames(spRPK))]
colnames(sp_F_mat) <- gsub("(FRNA_)*", "\\2", colnames(sp_F_mat))
FRNA.sizefactor <- SizeFactorCal(sp_F_mat)
sp_F_mat <- t(t(sp_F_mat) / FRNA.sizefactor)

sp_L_mat <- spRPK[,grepl("LRNA",colnames(spRPK))]
colnames(sp_L_mat) <- gsub("(LRNA_)*", "\\2", colnames(sp_L_mat))
LRNA.sizefactor <- SizeFactorCal(sp_L_mat)
sp_L_mat <- t(t(sp_L_mat) / LRNA.sizefactor)

saveRDS(FRNA.sizefactor, "../data/FRNA.sizefactor.RData")
saveRDS(LRNA.sizefactor, "../data/LRNA.sizefactor.RData")

## spike-ins table, for labelel rate estimation
SampleSpCounts <- NULL
for(i in seq_len(ncol(sp_F_mat))){
  SampleSpCounts <- rbind(SampleSpCounts,
                          data.frame(FRNA = sp_F_mat[,i],
                                     LRNA = sp_L_mat[,i],
                                     Sample = colnames(sp_L_mat)[i],
                                     SpikeIns = rownames(sp_F_mat),
                                     W = c(1, 0.1, 1, 0.1, 1),
                                     R = c(1, 1, 0.1, 0.1, 0) ))
}
saveRDS(SampleSpCounts, "../data/SampleSpikeCounts.RData")

# normalise tx tpm with spike-in RNA size factor
tx_F_mat <- txRPK[, grepl("FRNA", colnames(txRPK))]
colnames(tx_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tx_F_mat))
tx_F_mat <- sweep(tx_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tx_L_mat <- txRPK[, grepl("LRNA", colnames(txRPK))]
colnames(tx_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tx_L_mat))
tx_L_mat <- sweep(tx_L_mat, 2, LRNA.sizefactor, '/')

# write.table(tx_F_mat, "data/tx_FRNA_RPK_norm.txt",
#             quote = F, row.names = T, col.names = T, sep = "\t")
# write.table(tx_L_mat, "data/tx_LRNA_RPK_norm.txt",
#             quote = F, row.names = T, col.names = T, sep = "\t")
saveRDS(tx_L_mat, "../data/tx_L_mat.RData")
saveRDS(tx_F_mat, "../data/tx_F_mat.RData")

# normalise non-coding TU tpm with spike-in RNA size factor
tu_F_mat <- tuRPK[, grepl("FRNA", colnames(tuRPK))]
colnames(tu_F_mat) <- gsub("FRNA_(.*)","\\1", colnames(tu_F_mat))
tu_F_mat <- sweep(tu_F_mat, 2, FRNA.sizefactor, '/') # divide spike-in size factors

tu_L_mat <- tuRPK[, grepl("LRNA",colnames(tuRPK))]
colnames(tu_L_mat) <- gsub("LRNA_(.*)","\\1", colnames(tu_L_mat))
tu_L_mat <- sweep(tu_L_mat, 2, LRNA.sizefactor, '/')

saveRDS(tu_F_mat, "data/tu_FRNA_RPK_norm.rds")
saveRDS(tu_L_mat, "data/tu_LRNA_RPK_norm.rds")

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
