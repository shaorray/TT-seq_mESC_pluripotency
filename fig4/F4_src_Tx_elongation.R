# load data  -------------------------------------------------------------------------------
sample_conditions <- c("SL", "2i_2d", "mTORi_1d")
L_RNA_sf <- readRDS("../data/LRNA.sizefactor.RC.RData")
# gene.gr <- readRDS("../data/gene.gr.RData") 
# gene.gr <- gene.gr[width(gene.gr) > 2000 & width(gene.gr) < 1000000]

# use txdb gene references
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(gene.gr$gene_id),
                       keytype = "ENTREZID",
                       columns = "GENEID")
gene.gr$gene_id <- res$GENEID[match(gene.gr$gene_id, res$ENTREZID)]

# keep only active genes, longer than 2kb
gene.gr <- gene.gr[width(gene.gr) > 2000 & width(gene.gr) < 1000000]
gene.gr <- gene.gr[seqnames(gene.gr) %in% paste0("chr", c(1:19, "X", "Y"))]
txRPK <- readRDS("../data/txRPK_SL_2i.RData")
gene.gr <- sort(gene.gr[na.omit(match(rownames(txRPK)[rowMeans(txRPK) > 0.2], gene.gr$gene_id))])
names(gene.gr) <- gene.gr$gene_id # 10674 genes
# export.gff3(gene.gr, "../data/mm9.active.gene.gff3")

# plot comparison  -------------------------------------------------------------------------
simple_cor_plot <- function(x, y, viridis_type = "A", color_end = 1, .xlab, .ylab)
{
  data.frame(x = log(x), y = log(y)) %>%
    dplyr::filter(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)) %>%
    dplyr::filter(x > quantile(x, 0.05) & y > quantile(y, 0.05)) %>%
    dplyr::filter(x < quantile(x, 0.99) & y < quantile(y, 0.99)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_hex(bins= 50) +
    scale_fill_viridis_c(option = viridis_type, direction = -1, end = color_end) +
    xlab(.xlab) + ylab(.ylab) +
    theme_setting +
    theme(legend.position = "none")
}

# process coverages ------------------------------------------------------------------------
# ChIP coverage ----------------------------------------------------------------------------
# estimate MINUTE ChIP sample sizes
pol2s5p_input_size_factor = .countBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                                             pattern = "P1_IN_.*ALL.*bam$", full.names = T),
                                      gene.gr) %>% SizeFactorCal()
saveRDS(pol2s5p_input_size_factor, "data/pol2s5p_input_size_factor.RData")

# read Pol2s5p gene body + flanks (sandwich) coverage
bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                       pattern = "P1.*Pol.*ALL.*bam$", full.names = T)

pol2s5p_mat_sandwich <- readBam(bam_files = bam_files,
                                intervals = gene.gr,
                                pair_end = T,
                                stranded = F,
                                flanks = c(2000, 4000),
                                new_lens = c(20, 200, 40))

# normalise ChIP coverage between samples
pol2s5p_mat_sandwich_list <- list()
for (i in c(3, 1, 4)) {
  pol2s5p_mat_sandwich_list <- c(pol2s5p_mat_sandwich_list,
                                 list(pol2s5p_mat_sandwich[[i]] / pol2s5p_input_size_factor[i]))
}
names(pol2s5p_mat_sandwich_list) <- sample_conditions
# saveRDS(pol2s5p_mat_sandwich_list, "data/pol2s5p_gene_sandwich_mat.RData")


# TTseq coverage  ----------------------------------------------------------------------------
bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                       pattern = "LRNA.*(SL_|_2i_2d|mTORi_).*bam$", full.names = T)

TT_gene_sandwich_mat <- readBam(bam_files = bam_files,
                                intervals = gene.gr,
                                pair_end = T,
                                stranded = F, # include upstream reads on the main gene coverage
                                flanks = c(2000, 4000),
                                new_lens = c(20, 200, 40))

names(TT_gene_sandwich_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1", names(TT_gene_sandwich_mat))

# normalise TTseq coverage by spike-in size factors
for (i in names(TT_gene_sandwich_mat)) {
  TT_gene_sandwich_mat[[i]] <- TT_gene_sandwich_mat[[i]] / L_RNA_sf[i]
}
# averaging replicates
TT_gene_sandwich_mat_norm <- list()
for (i in c("SL", "2i", "mTORi")) {
  idx <- grep(i, names(TT_gene_sandwich_mat))
  TT_gene_sandwich_mat_norm <- c(TT_gene_sandwich_mat_norm,
                                 list(Reduce("+", TT_gene_sandwich_mat[idx]) / length(idx)))
}

names(TT_gene_sandwich_mat_norm) <- sample_conditions
# saveRDS(TT_gene_sandwich_mat_norm, "data/TT_gene_sandwich_mat_norm.RData")

TT_Pol2s5p_sandwich_mat <- list()
for (i in sample_conditions) {
  TT_Pol2s5p_sandwich_mat <- c(TT_Pol2s5p_sandwich_mat,
                          list(log1p(TT_gene_sandwich_mat_norm[[i]]) - log1p(pol2s5p_mat_sandwich_list[[i]])) )
}
names(TT_Pol2s5p_sandwich_mat) <- sample_conditions

# saveRDS(TT_Pol2s5p_sandwich_mat, "data/TT_Pol2s5p_sandwich_mat.RData")

# process reads density -------------------------------------------------------------------------------
# gene body density ------------------------------------------------------------------------------------------------------
# Pol2s5p gene body normed RPK
# gene.gr <- TU.DE.mm9.gr[TU.DE.mm9.gr$location == "protein_coding" & !is.na(TU.DE.mm9.gr$gene_id)]

gene_body_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                                  "P1.*Pol5p.*ALL.*bam$", full.names = T),
                           intervals = gene.gr,
                           stranded = F, paired.end = "ignore")
gene_body_mat <- gene_body_mat / width(gene.gr) * 1e3 # convert to reads density

pol2s5p_input_size_factor <- readRDS("data/pol2s5p_input_size_factor.RData")
gene_body_mat <- sweep(gene_body_mat, 2, pol2s5p_input_size_factor, "/")

gene_body_mat <- gene_body_mat[, c(3, 1, 4)]
rownames(gene_body_mat) <- gene.gr$gene_id
colnames(gene_body_mat) <- c("Pol2_SL", "Pol2_2i", "Pol2_mTORi")
gene_body_mat <- log10(gene_body_mat)
gene_body_mat <- gene_body_mat[!apply(gene_body_mat, 1,
                                      function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# TT-seq LRNA density  ----------------------------------------------------------------------------
bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                       pattern = "LRNA.*(SL_|_2i_2d|mTORi_2d).*(rep1|rep2).*bam$", full.names = T)

TT_gene_density_mat <- .countBam(bam_files = bam_files,
                                 intervals = gene.gr,
                                 stranded = F,
                                 paired.end = "ignore")

names(TT_gene_density_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1", names(TT_gene_density_mat))

# normalise TTseq coverage by spike-in size factors
TT_gene_density_mat <- sweep(TT_gene_density_mat, 2, L_RNA_sf[names(TT_gene_density_mat)], "/") 

# averaging replicates
tmp <- NULL
for (i in c("SL", "2i", "mTORi")) {
  idx <- grep(i, colnames(TT_gene_density_mat))
  tmp <- cbind(tmp,
               TT_gene_density_mat[, idx] %>%
                 "/"(width(gene.gr) / 1e3) %>% 
                 rowMeans() %>% log10() )
}
TT_gene_density_mat <- `rownames<-`(tmp, gene.gr$gene_id)
colnames(TT_gene_density_mat) <- c("Tx_RPK_SL", "Tx_RPK_2i", "Tx_RPK_mTORi")

# combine Pol2 and RNA synthesis
gene_body_production_mat <- cbind(gene_body_mat, TT_gene_density_mat[rownames(gene_body_mat), ])

gene_body_production_mat$LRNA_RPK <- (TU.DE.mm9.gr$baseMean_LRNA_sp / width(TU.DE.mm9.gr) * 1e3)[
  match(rownames(gene_body_production_mat), TU.DE.mm9.gr$gene_id)] %>% log10()

gene_body_production_mat <- gene_body_production_mat[
  !apply(gene_body_production_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# load 2014 Jonker et.al Fig3 data ------------------------------------------------------------------------
paper_est_speed <- read.table("../data/elife-02407-fig3-data1-v1", sep = "\t", header = 1)

gene_mtch <- biomaRt::getBM(filters = "refseq_mrna", attributes = c("refseq_mrna", "ensembl_gene_id"),
                      values = paper_est_speed$Genename, 
                      mart = biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl")))

paper_est_speed <- paper_est_speed[match(gene_mtch$refseq_mrna, paper_est_speed$Genename), ]
paper_est_speed$gene_id <- gene_mtch$ensembl_gene_id
paper_est_speed$speed <- paper_est_speed$Rate..bp.min. / 1e3

