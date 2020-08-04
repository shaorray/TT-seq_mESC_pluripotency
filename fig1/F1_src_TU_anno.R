setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# functions ------------------------------------------------------------------------
Count_RNA_reads <- function(TU.list, bam_files)
{
  require(Rsubread)
  locations <- unique(TU.list[[1]]$location)
  TU.out <- GRanges()
  for (loc in locations)
  {
    TU.tmp <- reduce(
      Reduce(c, lapply(TU.list, function(TU) TU[TU$location == loc]) )
    )
    TU.tmp$location <- loc
    TU.out %c=% TU.tmp
  }
  # append spike-in normalised read counts to TU GRanges output
  mcols(TU.out) %c=% featureCounts(bam_files,
                                   annot.ext = data.frame(GeneID = seq_along(TU.out),
                                                          Chr = seqnames(TU.out),
                                                          Start = start(TU.out),
                                                          End = end(TU.out),
                                                          Strand = strand(TU.out)),
                                   isPairedEnd=TRUE, nthreads = 8)$counts
  return(TU.out)
}

setNamesLRNACounts <- function(readCounts)
{
  readCounts <- readCounts[, !grepl("FRNA", colnames(readCounts))]
  colnames(readCounts) <- gsub("\\.Aligned.*", "\\1", colnames(readCounts))
  colnames(readCounts) <- gsub("\\.", "\\_", colnames(readCounts))
  colnames(readCounts) <- gsub("SLmTORi_24|SLmt_24", "mTORi_1d", colnames(readCounts))
  colnames(readCounts) <- gsub("mTORi_1d$", "mTORi_1d_rep2", colnames(readCounts))
  colnames(readCounts) <- gsub("SLmt48h", "mTORi_2d", colnames(readCounts))
  colnames(readCounts) <- gsub("LRNA\\_", "\\2", colnames(readCounts))
  # colnames(readCounts) <- gsub("^2i_", "2i_2d_", colnames(readCounts))
  # colnames(readCounts) <- gsub("2i7d_", "2i_7d_", colnames(readCounts))
  readCounts
  # readCounts <- readCounts[, names(LRNA.sizefactor)]
  # as.data.frame(t(t(as.matrix(readCounts)) / LRNA.sizefactor))
}

convert_expr <- function (TU.gr, TU.expr, bin_size = 200L, type = "TPM")
{ # normalise TU annotation expression in bins
  if (type == "TPM")
  {
    return(as.numeric(TU.expr) / sum(as.numeric(TU.expr)) * 1e6)
  }
  if (type == "RPKM")
  {
    return(as.numeric(TU.expr) / 1000 / sum(width(TU.gr) * as.numeric(TU.expr)) * 1e6)
  }
  if (type == "Read")
  {
    return(as.numeric(TU.expr) * width(TU.gr) / bin_size)
  }
}

jaccard_index_pair = function(TU.gr.list)
{
  locations <- unique(TU.list[[1]]$location)
  idx_mat_list <- list()
  for (l in locations) {
    n <- length(TU.gr.list)
    TU.tpm.list <- lapply(TU.gr.list, function(x) x[x$location == l])
    idx_mat <- expand.grid(1:n, 1:n)
    idx_mat <- idx_mat[!duplicated(t(apply(idx_mat, 1, sort))) & (idx_mat[,1] - idx_mat[,2]) != 0, ]
    idx_mat$Jaccard <- NA
    for (i in seq_len(nrow(idx_mat))) {
      idx_mat$Jaccard[i] <-
        sum(width(intersect(TU.tpm.list[[ idx_mat[i, 1] ]], TU.tpm.list[[ idx_mat[i, 2] ]]))) /
        sum(width(reduce(c(TU.tpm.list[[ idx_mat[i, 1] ]], TU.tpm.list[[ idx_mat[i, 2] ]]))))
    }
    idx_mat$Var1 <- factor(names(TU.gr.list)[idx_mat$Var1], levels = rev(names(TU.gr.list)))
    idx_mat$Var2 <- factor(names(TU.gr.list)[idx_mat$Var2], levels = names(TU.gr.list))
    idx_mat_list <- c(idx_mat_list, list(idx_mat))
  }
  names(idx_mat_list) <- locations
  idx_mat_list
}

# read in annotations ------------------------------------------------------------------------
TU.list <- list("SL_rep1" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SL_rep1.Aligned.sortedByCoord.out.gtf"),
                "SL_rep2" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SL_rep2.Aligned.sortedByCoord.out.gtf"),
                "SL_rep3" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SL_rep3.Aligned.sortedByCoord.out.gtf"),
                "2i_rep1" = importRanges("../data/TU_anno/other/TU_filter+LRNA_2i_rep1.Aligned.sortedByCoord.out.gtf"),
                "2i_rep2" = importRanges("../data/TU_anno/other/TU_filter+LRNA_2i_rep2.Aligned.sortedByCoord.out.gtf"),
                "mTORi_1d_rep1" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SLmTORi_1d_rep1.Aligned.sortedByCoord.out.gtf"),
                "mTORi_1d_rep2" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SLmTORi_1d_rep2.Aligned.sortedByCoord.out.gtf"),
                "mTORi_2d_rep1" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SLmTORi_2d_rep1.Aligned.sortedByCoord.out.gtf"),
                "mTORi_2d_rep2" = importRanges("../data/TU_anno/other/TU_filter+LRNA_SLmTORi_2d_rep2.Aligned.sortedByCoord.out.gtf"))

# count reads on annotated TUs
bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                        pattern = ".bam$", full.names = T)

# LRNA.sizefactor <- readRDS("../data/LRNA.sizefactor.RData")

TU.counts <- Count_RNA_reads(TU.list = TU.list, bam_files = bam_files)
mcols(TU.counts) <- cbind(mcols(TU.counts)[, 1],
                          setNamesLRNACounts(mcols(TU.counts)[, -1]))

# Liftover to mm9 for counting bam files
TU.list.mm9 <- lapply(TU.list, function(x)
  liftOver(x, chain = import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>%
    unlist)

TU.counts.mm9 <- Count_RNA_reads(
  TU.list = TU.list.mm9,
  bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                         pattern = ".bam$", full.names = T))
# saveRDS(TU.counts.mm9, "data/TU.read.counts.mm9.RData")

TU.TPM.mm9 <- TU.counts.mm9
mcols(TU.TPM.mm9)[, -1] <- apply(mcols(TU.counts.mm9)[, -1], 2, convert_expr, TU.gr = TU.counts.mm9)
mtch <- findOverlaps(TU.TPM.mm9, TU.list.mm9[[1]], ignore.strand = F)
TU.TPM.mm9$gene_id <- NA
TU.TPM.mm9$gene_id[queryHits(mtch)] <- TU.list.mm9[[1]]$gene_id[subjectHits(mtch)]

# save output
# saveRDS(TU.list, "data/TU.list.RData")
saveRDS(TU.counts, "data/TU.counts.RData")
saveRDS(TU.TPM.mm9, "data/TU.TPM.mm9.RData")

# export ncRNA sequences, for alignment free counting with all tx references
# will be applied for half-life calculation at the later steps
TU.mm9.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
                     TU.TPM.mm9)
names(TU.mm9.seq) = paste(TU.TPM.mm9$location,
                          seqnames(TU.TPM.mm9),
                          start(TU.TPM.mm9),
                          end(TU.TPM.mm9),
                          strand(TU.TPM.mm9), sep = "_")
TU.mm9.seq <- TU.mm9.seq[TU.TPM.mm9$location != "protein_coding"]
# writeXStringSet(TU.mm9.seq, "data/combined_ncRNA_seq.mm9.fa")

# reads distribution across TU types ------------------------------------------------------------------------
nascent_reads.list <- lapply(TU.list,
                     function(TU) aggregate(convert_expr(TU, TU$expr, bin_size = 200, type = "Read"),
                                            list(TU$location), sum))
nascent_reads <- data.frame()
for (i in seq_along(nascent_reads.list))
{
  nascent_reads.list[[i]][, 2] <- nascent_reads.list[[i]][, 2] / sum(nascent_reads.list[[i]][, 2]) * 100
  nascent_reads <- rbind(nascent_reads, cbind(nascent_reads.list[[i]], names(nascent_reads.list)[i]))
}
colnames(nascent_reads) <- c("location", "Abundance", "Sample")

# keep highly consistent TU types
tu_of_interest <- c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA")
nascent_reads$location[nascent_reads$location == "protein_coding"] <- "mRNA"
nascent_reads <- nascent_reads[nascent_reads$location %in% tu_of_interest, ]
nascent_reads$location <- factor(nascent_reads$location, levels = tu_of_interest)
nascent_reads$Abundance <- as.numeric(nascent_reads$Abundance)

write.table(nascent_reads, "data/Fig1_Tx_read_abundance.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

# total reads, FRNA
total_reads <- TU.counts %>% mcols() %>% as.data.frame() %>%
  dplyr::select(matches("FRNA")) %>%
  dplyr::group_by(TU.counts$location) %>%
  summarise_all(sum) %>%
  as.data.frame() # sum reads for each type
total_reads[, -1] <- total_reads[, -1] %>%
  sapply(function(x) x / sum(x) * 100 ) %>% unlist() # convert to percent
total_reads[, 1][total_reads[, 1] == "protein_coding"] <- "mRNA"
total_reads <- total_reads[total_reads[, 1] %in% tu_of_interest, ]
total_reads <- reshape2::melt(total_reads,
                              id = "TU.counts$location")
colnames(total_reads) <- c("location", "Sample", "Abundance")
total_reads$Sample <- gsub("\\.", "_", gsub("FRNA.(.*).Aligned.*", "\\1", total_reads$Sample))
total_reads$location <- factor(total_reads$location, levels = tu_of_interest)
total_reads$Abundance <- as.numeric(total_reads$Abundance)

write.table(total_reads, "data/Fig1_Total_read_abundance.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")


# TPM ------------------------------------------------------------------------
total_TPM <- data.frame()
for (i in seq_along(TU.list))
{
  total_TPM <- rbind(total_TPM,
                    cbind(TU.list[[i]]$location,
                          convert_expr(TU.list[[i]]$expr, bin_size = 200, type = "TPM"),
                          names(TU.list)[i])
                    )
}
colnames(total_TPM) <- c("location", "TPM", "Sample")
total_TPM$location <- as.character(total_TPM$location)
total_TPM$location[total_TPM$location == "protein_coding"] <- "mRNA"
total_TPM <- total_TPM[total_TPM$location %in%
                           c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA"), ]
total_TPM$location <- factor(total_TPM$location ,
                              levels = c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA"))
total_TPM$TPM <- as.numeric(as.character(total_TPM$TPM))
saveRDS(total_TPM, "data/total_TPM.RData")

# ncRNA DE analysis ------------------------------------------------------------------------
pacman::p_load(DESeq2)
# mm10 TU read counts
count_table <- mcols(TU.counts)[, -1]
dds <- DESeqDataSetFromMatrix(round(as.matrix(count_table)),
                              colData = data.frame(condition = gsub("_rep.", "\\2", colnames(mcols(TU.counts))[-1]) ),
                              design = ~ condition)
dds <- DESeq(dds)
res_2i <- results(dds, contrast = c("condition", "2i_2d", "SL"))
res_mTORi <- results(dds, contrast = c("condition", "mTORi_1d", "SL"))

TU.DE <- TU.counts
mcols(TU.DE) <- data.frame("location" = mcols(TU.DE)[, 1],
                           "baseMean" = res_2i[, "baseMean"],
                           "log2FoldChange_2i" = res_2i[, "log2FoldChange"],
                           "padj_2i" = res_2i[, "padj"],
                           "log2FoldChange_mTORi" = res_mTORi[, "log2FoldChange"],
                           "padj_mTORi" = res_mTORi[, "padj"]
                           )
# saveRDS(TU.DE, "data/TU.DE.RData")

# add attributes, direction, enhancer ----------------------------------------------------------------------------
if (T) {
  # enhancer in development
  F5_enhancer <- importRanges("../data/F5.mm10.enhancers.bed")
  
  TU.DE.intergenic <- TU.DE[TU.DE$location == "intergenic" & TU.DE$baseMean > 0]
  TU.intergenic_as <- promoters(TU.DE.intergenic, upstream = 500, downstream = 0)
  levels(strand(TU.intergenic_as)) <- c("-", "+", "*")
  
  mtch_dir <- findOverlaps(TU.DE.intergenic, TU.intergenic_as)
  TU.DE.intergenic$direction <- ifelse(mtch_dir %>% countQueryHits() > 0,
                                       "Bidirectional", "Unidirectional" )
  TU.DE.intergenic$enhancer <- ifelse(findOverlaps(F5_enhancer, TU.DE.intergenic) %>% countSubjectHits > 0,
                                      "Enhancer", "TX")
  TU.DE.intergenic$enhancer_direction <- paste(TU.DE.intergenic$direction, TU.DE.intergenic$enhancer, sep = "_")
  # major / minor direction by expr
  mtch_dir <- as.matrix(mtch_dir)
  for( i in seq_len(nrow(mtch_dir)) )
  { # sort by reads density
    if (which.max(TU.DE.intergenic$baseMean[mtch_dir[i, ]] / width(TU.DE.intergenic[mtch_dir[i, ]])) == 2)
      mtch_dir[i, ] <- rev(mtch_dir[i, ])
  }
  mtch_dir <- mtch_dir[!duplicated(mtch_dir), ]
  
  TU.DE.intergenic$main_direction <- TRUE
  TU.DE.intergenic$main_direction[mtch_dir[, 2]] <- FALSE
  
  # chromHMM
  ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
  ChromHMM$name <- gsub(".*_", "\\2", ChromHMM$name)
  TU.DE.intergenic$chromHMM <- ChromHMM$name[findOverlaps(promoters(TU.DE.intergenic, upstream = 1, downstream = 0),
                                                          ChromHMM) %>% subjectHits]
  
  # TE
  rmsk <- importRanges("/mnt/0E471D453D8EE463/genomeDir/rmsk/mm10_rmsk_TE.gtf")
  mtch_TE <- findOverlaps(TU.DE.intergenic, rmsk)
  TU.DE.intergenic$TE_class <- "None"
  TU.DE.intergenic$TE_family <- "None"
  TU.DE.intergenic$TE_class[queryHits(mtch_TE)] <- rmsk$class_id[subjectHits(mtch_TE)]
  TU.DE.intergenic$TE_family[queryHits(mtch_TE)] <- rmsk$family_id[subjectHits(mtch_TE)]
  
  # evo conservation
  phastCons.list <- list.files('/mnt/0E471D453D8EE463/genomeDir/phastCons', full.names = T)
  phastConsnMat <- sapply(phastCons.list,
                          function (bw_file)
                            summary(BigWigFile(bw_file),
                                    TU.DE.intergenic,
                                    defaultValue = 0,
                                    as = "matrix") )
  colnames(phastConsnMat) <- gsub(".*(mm10.*).bw", "\\1", colnames(phastConsnMat))
  mcols(TU.DE.intergenic) <- cbind(mcols(TU.DE.intergenic), phastConsnMat)
  
  # ChromHMM state coeffiencts on log2FC, uni- / bi-direction intergenic TU
  intergenic_2i_DE <- mcols(TU.DE.intergenic)[TU.DE.intergenic$enhancer == "TX", -c(1,2,5,6,9)] %>% as.data.frame
  coef_dat <- rbind(
    cbind("Bidirection", summary(lm( log2FoldChange_2i ~ chromHMM,
                                     intergenic_2i_DE[intergenic_2i_DE$direction == "Bidirectional", ]))$coef),
    cbind("Unidirection", summary(lm( log2FoldChange_2i ~ chromHMM,
                                      intergenic_2i_DE[intergenic_2i_DE$direction == "Unidirectional", ]))$coef)
  ) %>% as.data.frame
  coef_dat$chromHMM <- gsub(".1", "\\2", gsub("(chromHMM|X.)", "\\2", rownames(coef_dat)))
  coef_dat$chromHMM[grep("Intercept", coef_dat$chromHMM)] <- "ActivePromoter"
  coef_dat$chromHMM <- factor(coef_dat$chromHMM, levels = sort(unique(coef_dat$chromHMM), T) )
  coef_dat <- coef_dat[!grepl("Transcription", coef_dat$chromHMM), ]
  coef_dat$Estimate <- as.numeric(as.character(coef_dat$Estimate))
  coef_dat$"Std. Error" <- as.numeric(as.character(coef_dat$"Std. Error"))
  saveRDS(coef_dat, "../fig3/data/coef_dat.RData")
  
  # Bidirectional intergenic TU
  bi_intergenic_DE <- with(TU.DE.intergenic[c(mtch_dir)],
                           data.frame(log2FoldChange = c(log2FoldChange_2i, log2FoldChange_mTORi),
                                      padj = c(padj_2i, padj_mTORi),
                                      Sample = c(rep("2i", length(mtch_dir)),
                                                 rep("mTORi", length(mtch_dir) ) ),
                                      Direction = c(rep("Major", nrow(mtch_dir)),
                                                    rep("Minor", nrow(mtch_dir))),
                                      Enhancer = rep(enhancer, 2) ))
  saveRDS(bi_intergenic_DE, "../fig3/data/bi_intergenic_DE.RData")
  
  # uni directional intergenic TU DE
  uni_intergenic_DE <- with(TU.DE.intergenic[-c(mtch_dir)],
                            data.frame(log2FoldChange = c(log2FoldChange_2i, log2FoldChange_mTORi),
                                       padj = c(padj_2i, padj_mTORi),
                                       Sample = c(rep("2i", length(TU.DE.intergenic) - length(mtch_dir)),
                                                  rep("mTORi", length(TU.DE.intergenic) - length(mtch_dir))),
                                       Enhancer = rep(enhancer, 2) ))
  saveRDS(uni_intergenic_DE, "../fig3/data/uni_intergenic_DE.RData")
}

# mm9, spike-in norm DE ---------------------------------------------------------------------------
count_table <- mcols(TU.counts.mm9)[, -1]
count_table <- count_table[, !grepl("SL2i", colnames(count_table))]
L_sample_idx <- grepl("LRNA", colnames(count_table))
dds <- DESeqDataSetFromMatrix(round(as.matrix(count_table[, L_sample_idx])),
                              colData = data.frame(condition = gsub(".rep.*", "\\2", colnames(count_table)[L_sample_idx]) ),
                              design = ~ condition)
sizeFactors(dds) <- readRDS("../data/LRNA.sizefactor.RData")
dds <- DESeq(dds)
L_res_2i <- results(dds, contrast = c("condition", "LRNA.2i.2d", "LRNA.SL"))
L_res_mTORi <- results(dds, contrast = c("condition", "LRNA.mTORi.1d", "LRNA.SL"))

dds <- DESeqDataSetFromMatrix(round(as.matrix(count_table[, !L_sample_idx])),
                              colData = data.frame(condition = gsub(".rep.*", "\\2", colnames(count_table)[!L_sample_idx]) ),
                              design = ~ condition)
sizeFactors(dds) <- readRDS("../data/FRNA.sizefactor.RData")[1:10]
dds <- DESeq(dds)
F_res_2i <- results(dds, contrast = c("condition", "FRNA.2i.2d", "FRNA.SL"))
F_res_mTORi <- results(dds, contrast = c("condition", "FRNA.mTORi.1d", "FRNA.SL"))

TU.DE.mm9 <- TU.counts.mm9 # save as GRanges
mcols(TU.DE.mm9) <- data.frame("location" = mcols(TU.counts.mm9)[, 1],
                               "L_baseMean" = L_res_2i[, "baseMean"],
                               "L_log2FC_2i" = L_res_2i[, "log2FoldChange"],
                               "L_padj_2i" = L_res_2i[, "padj"],
                               "L_log2FC_mTORi" = L_res_mTORi[, "log2FoldChange"],
                               "L_padj_mTORi" = L_res_mTORi[, "padj"],
                               "F_baseMean" = F_res_2i[, "baseMean"],
                               "F_log2FC_2i" = F_res_2i[, "log2FoldChange"],
                               "F_padj_2i" = F_res_2i[, "padj"],
                               "F_log2FC_mTORi" = F_res_mTORi[, "log2FoldChange"],
                               "F_padj_mTORi" = F_res_mTORi[, "padj"] )

saveRDS(TU.DE.mm9, "data/TU.DE.external.norm.RData")

# DE with kallisto TPM --------------------------------------------------------------------------
txRPK <- readRDS("../data/txRPK_SL_2i.RData") # raw tmp without normalization
L_sample_idx <- grepl("LRNA", colnames(txRPK))
# LRNA
dds <- DESeqDataSetFromMatrix(round(as.matrix(txRPK[, L_sample_idx] * 100)), # make the values friendly to the nb model
                              colData = data.frame(condition = gsub(".rep.*", "\\2", colnames(txRPK)[L_sample_idx]) ),
                              design = ~ condition) 
sizeFactors(dds) <- readRDS("../data/LRNA.sizefactor.RData")
dds <- DESeq(dds)
L_res_2i <- results(dds, contrast = c("condition", "LRNA_2i_2d", "LRNA_SL"))
L_res_mTORi <- results(dds, contrast = c("condition", "LRNA_mTORi_1d", "LRNA_SL"))
# FRNA
dds <- DESeqDataSetFromMatrix(round(as.matrix(txRPK[, !L_sample_idx] * 100)),
                              colData = data.frame(condition = gsub(".rep.*", "\\2", colnames(txRPK)[!L_sample_idx]) ),
                              design = ~ condition) 
sizeFactors(dds) <- readRDS("../data/FRNA.sizefactor.RData")
dds <- DESeq(dds)
F_res_2i <- results(dds, contrast = c("condition", "FRNA_2i_2d", "FRNA_SL"))
F_res_mTORi <- results(dds, contrast = c("condition", "FRNA_mTORi_1d", "FRNA_SL"))

TX.DE.tpm <- data.frame("gene_id" = rownames(F_res_2i),
                        "L_baseMean" = L_res_2i[, "baseMean"],
                        "L_log2FC_2i" = L_res_2i[, "log2FoldChange"],
                        "L_padj_2i" = L_res_2i[, "padj"],
                        "L_log2FC_mTORi" = L_res_mTORi[, "log2FoldChange"],
                        "L_padj_mTORi" = L_res_mTORi[, "padj"],
                        "F_baseMean" = F_res_2i[, "baseMean"],
                        "F_log2FC_2i" = F_res_2i[, "log2FoldChange"],
                        "F_padj_2i" = F_res_2i[, "padj"],
                        "F_log2FC_mTORi" = F_res_mTORi[, "log2FoldChange"],
                        "F_padj_mTORi" = F_res_mTORi[, "padj"] )

saveRDS(TX.DE.tpm, "data/TX.DE.tpm.external.norm.RData")

