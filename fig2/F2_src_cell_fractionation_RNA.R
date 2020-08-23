# Rui 2020 Jan

# analyze cell fractionation RNA-seq dataset from https://www.nature.com/articles/nature20149

# clean read counts ---------------------------------------------------------------
file_paths <- list.files('../../data/Fractionation_output', pattern = 'abundance.tsv', full.names = T, recursive = T)
sample_names <- sapply(file_paths, function(a) unlist(strsplit(a, '/'))[5] ) %>% unname

# read kallisto counts
results <- SummarizedExperiment::readKallisto(file_paths, as = 'matrix')
colnames(results) <- sample_names

# remove inactive genes
results <- results[apply(results, 1, function(x) any(x > 0) ), ]

# keep total RNA non-zero, remove alternative txs of ptrotein coding genes
results <- results[apply(results[, c(2,4,6)], 1, function(x) any(x > 0) & mean(x) > 1), ]
results <- keepOneTx(count_table = results, rowname_gene_id = T)

F_norm_results <- t(t(results) / SizeFactorCal(results)) %>% log1p
# F_norm_results <- F_norm_results / rowMax(F_norm_results)
colnames(F_norm_results) <- gsub("2016_Jesse_", "\\1", colnames(F_norm_results))

saveRDS(F_norm_results, 'data/Fractionation_RNA_log_normed.RData')
write.table(F_norm_results, 'data/Fractionation_RNA_log_normed.txt',
            sep = '\t', quote = F, row.names = T)

# clustering test ------------------------------------------------------------------------
p = prcomp(F_norm_results, center = T, scale. = T)
# smoothScatter(p$rotation[, 1:2])
clusters <- kmeans( t(scale(t(F_norm_results))), 5, nstart = 20 )
par(mfrow = c(3,2))
for( i in 1:5 )
{
  barplot(colMeans(t(scale(t(F_norm_results[clusters$cluster==i, ])))), las = 2, col = colors_9[i])
  title(i, adj=0, line = -0.3)
}
plot(p$x[, 1:2], pch = 19, cex = 0.2, col = colors_9[clusters$cluster])

# enrichment annotations ---------------------------------------------------------------
# evaluate RNA poly adenylylation, subcellular enrichment and expression levels
RNA_types <- c('PolyA', 'Non_polyA')
locations <- c('Chromatin', 'Cytoplasm', 'Nucleoplasm')
expression_states <- c('active', 'weak')

gene_attributes <- data.frame(gene_id = rownames(F_norm_results))

gene_attributes$RNA_type <- RNA_types[1]
gene_attributes$RNA_type[rowMeans(F_norm_results[, c(2,4,6)]) > rowMeans(F_norm_results[, -c(2,4,6)])] <- RNA_types[2]
gene_attributes$RNA_polyA_enrichment <- rowMeans(F_norm_results[, -c(2,4,6)]) - rowMeans(F_norm_results[, c(2,4,6)])

location_table <- cbind(rowMeans(F_norm_results[, 1:2]), rowMeans(F_norm_results[, 3:4]), rowMeans(F_norm_results[, 5:6]))
location_assigned <- apply(location_table, 1, which.max)
gene_attributes$subcellular_location <- locations[location_assigned]
gene_attributes$subcellular_enrichment <- location_table[matrix(c(seq_len(nrow(gene_attributes)), location_assigned), ncol = 2)] / rowMeans(location_table)

gene_attributes$expression_state <- expression_states[2 - (rowMeans(F_norm_results) > 3)]
gene_attributes$mean_expression <- rowMeans(F_norm_results)

# write.table(gene_attributes, 'data/Fractionation_RNA_gene_attributes.txt', sep = '\t', quote = F, row.names = F)
saveRDS(gene_attributes, 'data/Fractionation_gene_attributes.RData')

# preprocessing RNA labeled/decay rates ------------------------------------------------------
tx_label_rate <- readRDS('../figS2/data/sample_Tx_counts_Rates_combined.RData')
tx_label_rate <- tx_label_rate[tx_label_rate$Sample == 'SL', ]
tx_label_rate <- tx_label_rate[which(tx_label_rate$FRNA > 0 &
                                 tx_label_rate$LRNA > 0 & 
                                 tx_label_rate$Labeled_rate > 0 & 
                                 tx_label_rate$Labeled_rate < 2), ] %>% as.data.frame
saveRDS(tx_label_rate, "data/SL_tx_label_rate.RData")
# filter fractionation data
F_norm_results <- F_norm_results[rowSums(F_norm_results) > 10, ] # keep active genes
F_norm_results <- F_norm_results / rowSums(F_norm_results) # normalise gene expression

# make data collection
fraction_label_rate_dat <- F_norm_results[ rownames(F_norm_results) %in% tx_label_rate$gene_id, ]
fraction_label_rate_dat <- cbind(fraction_label_rate_dat,
                                 "r" = tx_label_rate$Labeled_rate[match(rownames(fraction_label_rate_dat), tx_label_rate$gene_id)]
) %>% as.data.frame() # 10932  genes

# read in RNA decay rate estimation from SLAM-seq, Herzog et.al, 2017
pacman::p_load(gdata)
SLAM.table = read.xls(xls = "../data/41592_2017_BFnmeth4435_MOESM4_ESM.xls", sheet = 1)

# convert gene id
ks <- keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "ENSEMBL")
res <- biomaRt::select(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL",
                       columns = c("ENTREZID", "SYMBOL"))
SLAM.table$gene_id <- res$ENSEMBL[ match(SLAM.table$Name, res$SYMBOL) ]
SLAM.table <- SLAM.table[!is.na(SLAM.table$gene_id) & SLAM.table$k..cpm.h. > 0, ] # 6534 genes

# combine existing gene counts as the explantory virables
same_genes <- intersect.Vector(intersect.Vector(SLAM.table$gene_id,
                                                rownames(F_norm_results)),
                               as.character(tx_label_rate$gene_id)) # 5018

fit_data <- cbind(F_norm_results[match(same_genes, rownames(F_norm_results)), ], 
                  "k" = SLAM.table$k..cpm.h.[match(same_genes, SLAM.table$gene_id)] / 60,
                  "r" = tx_label_rate$Labeled_rate[match(same_genes, tx_label_rate$gene_id)]) %>% as.data.frame()
# convert 5 minutes labeled rates to theoretical decay rates, with the equation Y(t) = α / β * (1 - e^(-t*β))
fit_data$r <- log(1 - fit_data$r) / (-5) 
fit_data <- fit_data[!apply(fit_data, 1, function(x) any(is.na(x))), ] # 4993 genes

saveRDS(fit_data, "data/fractionation_fit_data.RData")
