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

setNames_normCounts <- function(readCounts, LRNA.sizefactor, FRNA.sizefactor)
{
  # readCounts <- readCounts[, !grepl("FRNA", colnames(readCounts))]
  colnames(readCounts) <- gsub("\\.Aligned.*", "\\1", colnames(readCounts))
  colnames(readCounts) <- gsub("\\.", "\\_", colnames(readCounts))
  colnames(readCounts) <- gsub("SLmTORi_24|SLmt_24", "mTORi_1d", colnames(readCounts))
  colnames(readCounts) <- gsub("mTORi_1d$", "mTORi_1d_rep2", colnames(readCounts))
  colnames(readCounts) <- gsub("SLmt48h", "mTORi_2d", colnames(readCounts))
  # colnames(readCounts) <- gsub("LRNA\\_", "\\2", colnames(readCounts))
  names(FRNA.sizefactor) <- paste0("FRNA_", names(FRNA.sizefactor))
  names(LRNA.sizefactor) <- paste0("LRNA_", names(LRNA.sizefactor))
  
  round(cbind(t(t(as.matrix(readCounts[, names(FRNA.sizefactor)])) / FRNA.sizefactor),
              t(t(as.matrix(readCounts[, names(LRNA.sizefactor)])) / LRNA.sizefactor)), 2)
}

find_mm10_gene_id <- function(intervals) {
  gene.gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
  gene.gr <- gene.gr[gene.gr$gene_biotype == 'protein_coding']
  gene.gr <-  gene.gr[width(gene.gr) < 2500000]
  
  mtch <- findOverlaps(intervals, gene.gr, ignore.strand = F)
  gene_ids <- rep(NA, length(intervals))
  gene_ids[queryHits(mtch)] <- gene.gr$gene_id[subjectHits(mtch)]
  gene_ids
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
# mm10 annotations
TU.list <- list("SL" = importRanges("../data/TU_anno/mm10/TU_filter+LRNA_SL.gff3"),
                "2i_2d" = importRanges("../data/TU_anno/mm10/TU_filter+LRNA_2i_2d.gff3"),
                "mTORi_1d" = importRanges("../data/TU_anno/mm10/TU_filter+LRNA_mTORi_1d.gff3"),
                "mTORi_2d" = importRanges("../data/TU_anno/mm10/TU_filter+LRNA_mTORi_2d.gff3") )

# count reads on annotated TUs
bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                        pattern = ".bam$", full.names = T)

LRNA.sizefactor <- readRDS("../data/LRNA.sizefactor.RC.RData")
FRNA.sizefactor <- readRDS("../data/FRNA.sizefactor.RC.RData")

# featureCounts
TU.counts <- Count_RNA_reads(TU.list = TU.list, bam_files = bam_files)
mcols(TU.counts) <- cbind(mcols(TU.counts)[, 1],
                          setNames_normCounts(mcols(TU.counts)[, -1], LRNA.sizefactor, FRNA.sizefactor))

library(magrittr)
TU.counts.mat <- data.frame("seqname" = as.character(seqnames(TU.counts)), 
                            "start" = start(TU.counts),
                            "end" = end(TU.counts),
                            "strand" = as.character(strand(TU.counts)),
                            "location" = as.character(mcols(TU.counts)[, 1]),
                            "gene_id" = find_mm10_gene_id(TU.counts))
TU.counts.mat[TU.counts.mat$location != "protein_coding" & !is.na(TU.counts.mat$gene_id), "gene_id"] <- NA
TU.counts.mat <- cbind(TU.counts.mat, 
                       (mcols(TU.counts)[, -1]) %>% as.data.frame() %<>%
                         mutate_if(is.factor, as.character) %<>% 
                         mutate_if(is.character, as.numeric))

if (T) {
  # mm9 read counts
  TU.list.mm9 <- list("SL" = importRanges("../data/TU_anno/mm9/TU_filter+LRNA_SL.mm9.gff3"),
                      "2i_2d" = importRanges("../data/TU_anno/mm9/TU_filter+LRNA_2i_2d.mm9.gff3"),
                      "mTORi_1d" = importRanges("../data/TU_anno/mm9/TU_filter+LRNA_mTORi_1d.mm9.gff3"),
                      "mTORi_2d" = importRanges("../data/TU_anno/mm9/TU_filter+LRNA_mTORi_2d.mm9.gff3") )
  bam_files.mm9 <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                              pattern = ".bam$", full.names = T)
  TU.counts.mm9 <- Count_RNA_reads(TU.list.mm9, bam_files.mm9)

  FRNA.sizefactor <- readRDS("../fig1/data/sizefactor.list.mm9.RData")[["FRNA.sizefactor"]]
  LRNA.sizefactor <- readRDS("../fig1/data/sizefactor.list.mm9.RData")[["LRNA.sizefactor"]]
  
  mtch <- findOverlaps(TU.counts.mm9, TU.list.mm9[[1]], ignore.strand = F)
  TU.counts.mm9$gene_id <- NA
  TU.counts.mm9$gene_id[queryHits(mtch)] <- TU.list.mm9[[1]]$gene_id[subjectHits(mtch)]
  
  TU.counts.mat.mm9 <- data.frame("seqname" = as.character(seqnames(TU.counts.mm9)), 
                              "start" = start(TU.counts.mm9),
                              "end" = end(TU.counts.mm9),
                              "strand" = as.character(strand(TU.counts.mm9)),
                              "location" = mcols(TU.counts.mm9)[, 1],
                              "gene_id" = TU.counts.mm9$gene_id)
  TU.counts.mat.mm9[TU.counts.mat.mm9$location != "protein_coding" & !is.na(TU.counts.mat.mm9$gene_id), "gene_id"] <- NA
  TU.counts.mat.mm9 <- cbind(TU.counts.mat.mm9, 
                             setNames_normCounts(mcols(TU.counts.mm9)[, -1], LRNA.sizefactor, FRNA.sizefactor))
  mcols(TU.counts.mm9) <- cbind(mcols(TU.counts.mm9)[, 1],
                                setNames_normCounts(mcols(TU.counts.mm9)[, -1], LRNA.sizefactor, FRNA.sizefactor))

  
}

# export ncRNA sequences, for alignment free counting with all tx references
# for half-life and copy estimation in figure 2
if (F) {
  TU.mm10.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                        TU.counts)
  names(TU.mm10.seq) = paste(TU.counts$location,
                             seqnames(TU.counts),
                             start(TU.counts),
                             end(TU.counts),
                             strand(TU.counts), 
                             sep = "_")
  TU.mm10.seq <- TU.mm10.seq[TU.counts$location != "protein_coding"]
  writeXStringSet(TU.mm10.seq, "data/combined_ncRNA_seq.mm10.fa")
  
  seq_lengths <- BSgenome::getBSgenome("mm9") %>% 
    seqlengths() %>% '['(names(seqlengths(TU.counts.mm9)))
  
  for (i in names(seq_lengths)) {
    start(TU.counts.mm9[seqnames(TU.counts.mm9) == i & start(TU.counts.mm9) > seq_lengths[i]]) <- 
      seq_lengths[i]
    end(TU.counts.mm9[seqnames(TU.counts.mm9) == i & end(TU.counts.mm9) > seq_lengths[i]]) <-
      seq_lengths[i]
  }
  
  TU.mm9.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
                       TU.counts.mm9)
  names(TU.mm9.seq) = paste(TU.counts.mm9$location,
                            seqnames(TU.counts.mm9),
                            start(TU.counts.mm9),
                            end(TU.counts.mm9),
                            strand(TU.counts.mm9), 
                            sep = "_")
  TU.mm9.seq <- TU.mm9.seq[TU.counts.mm9$location != "protein_coding"]
  writeXStringSet(TU.mm9.seq, "data/combined_ncRNA_seq.mm9.fa")
}

# ----------------------------------------------------------------------------------------------------------
if (T) {
  TU.list <- list("SL_rep1" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_SL_rep1.Aligned.sortedByCoord.out.gff3"),
                  "SL_rep2" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_SL_rep2.Aligned.sortedByCoord.out.gff3"),
                  "SL_rep3" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_SL_rep3.Aligned.sortedByCoord.out.gff3"),
                  "2i_rep1" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_2i_2d_rep1.Aligned.sortedByCoord.out.gff3"),
                  "2i_rep2" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_2i_2d_rep2.Aligned.sortedByCoord.out.gff3"),
                  "mTORi_1d_rep1" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_mTORi_1d_rep1.Aligned.sortedByCoord.out.gff3"),
                  "mTORi_1d_rep2" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_mTORi_1d_rep2.Aligned.sortedByCoord.out.gff3"),
                  "mTORi_2d_rep1" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_mTORi_2d_rep1.Aligned.sortedByCoord.out.gff3"),
                  "mTORi_2d_rep2" = importRanges("../data/TU_anno/mm10_rep/TU_filter+LRNA_mTORi_2d_rep2.Aligned.sortedByCoord.out.gff3"))
}
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

saveRDS(nascent_reads, "data/nascent_reads.RData")

# TPM per location -------------------------------------------------------------------
if (T) {
  TU.counts2 <- mcols(TU.counts)[, -(1:11)] %>% as.data.frame() %<>%
    mutate_if(is.factor, as.character) %<>% 
    mutate_if(is.character, as.numeric) %>% "/"(width(TU.counts) * 1e3)
  TU.counts2 <- sweep(TU.counts2, 2, colSums(TU.counts2) / 1e6, "/")
  TU.counts2 <- cbind("location" = mcols(TU.counts)[, 1], 
                      as.data.frame(TU.counts2))
  
  total_TPM <- reshape::melt(TU.counts2)
  colnames(total_TPM) <- c("location", "Sample", "TPM")
  total_TPM$Sample <- gsub("X", "", total_TPM$Sample)
  total_TPM$location <- as.character(total_TPM$location)
  total_TPM$location[total_TPM$location == "protein_coding"] <- "mRNA"
  total_TPM$location[total_TPM$location == "antisense"] <- "asRNA"
  total_TPM <- total_TPM[total_TPM$location %in%
                           c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA"), ]
  total_TPM$location <- factor(total_TPM$location ,
                               levels = c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA"))
  total_TPM$TPM <- as.numeric(as.character(total_TPM$TPM))
  total_TPM$Sample <- factor(total_TPM$Sample, levels = unique(total_TPM$Sample)[c(7:9, 1:6)])
  saveRDS(total_TPM, "data/total_TPM.RData")
}

# DE analysis ----------------------------------------------------------------------------
library(DESeq2)
# mm10 TU LRNA read counts
# internal normalized or spike-in normalized differential expression
get_DE_res <- function(TU.counts.mat, TU.gr, .include_SL2i = F) {
  count_table_FRNA <- TU.counts.mat %>% dplyr::select(matches("FRNA"))
  count_table_LRNA <- TU.counts.mat %>% dplyr::select(matches("LRNA"))
  dds_FRNA <- DESeqDataSetFromMatrix(round(as.matrix(count_table_FRNA)),
                                     colData = data.frame(condition = gsub("FRNA_(.*)_rep.", "\\1", colnames(count_table_FRNA)) ),
                                     design = ~ condition)
  dds_LRNA <- DESeqDataSetFromMatrix(round(as.matrix(count_table_LRNA)),
                                     colData = data.frame(condition = gsub("LRNA_(.*)_rep.", "\\1", colnames(count_table_LRNA)) ),
                                     design = ~ condition)
  dds_FRNA <- DESeq(dds_FRNA)
  dds_LRNA <- DESeq(dds_LRNA)
  
  res_FRNA_2i <- results(dds_FRNA, contrast = c("condition", "2i_2d", "SL"))
  res_FRNA_mTORi <- results(dds_FRNA, contrast = c("condition", "mTORi_1d", "SL"))
  res_LRNA_2i <- results(dds_LRNA, contrast = c("condition", "2i_2d", "SL"))
  res_LRNA_mTORi <- results(dds_LRNA, contrast = c("condition", "mTORi_1d", "SL"))
  
  if (.include_SL2i) {
    res_FRNA_SL2i <- results(dds_FRNA, contrast = c("condition", "SL2i_2d", "SL"))
    res_LRNA_SL2i <- results(dds_LRNA, contrast = c("condition", "SL2i_2d", "SL"))
  }
  
  sizeFactors(dds_FRNA) <- rep(1, ncol(count_table_FRNA))
  sizeFactors(dds_LRNA) <- rep(1, ncol(count_table_LRNA))
  
  dds_FRNA <- DESeq(dds_FRNA)
  dds_LRNA <- DESeq(dds_LRNA)
  
  res_FRNA_2i_sp <- results(dds_FRNA, contrast = c("condition", "2i_2d", "SL"))
  res_FRNA_mTORi_sp <- results(dds_FRNA, contrast = c("condition", "mTORi_1d", "SL"))
  res_LRNA_2i_sp <- results(dds_LRNA, contrast = c("condition", "2i_2d", "SL"))
  res_LRNA_mTORi_sp <- results(dds_LRNA, contrast = c("condition", "mTORi_1d", "SL"))
  
  if (.include_SL2i) {
    res_FRNA_SL2i_sp <- results(dds_FRNA, contrast = c("condition", "SL2i_2d", "SL"))
    res_LRNA_SL2i_sp <- results(dds_LRNA, contrast = c("condition", "SL2i_2d", "SL"))
  }
  
  TU.DE <- TU.gr
  res_dat <- data.frame("location" = mcols(TU.DE)[, 1],
                        "gene_id" = TU.counts.mat$gene_id,
                        "baseMean_FRNA" = res_FRNA_2i[, "baseMean"],
                        "log2FoldChange_FRNA_2i" = res_FRNA_2i[, "log2FoldChange"],
                        "padj_FRNA_2i" = res_FRNA_2i[, "padj"],
                        "log2FoldChange_FRNA_mTORi" = res_FRNA_mTORi[, "log2FoldChange"],
                        "padj_FRNA_mTORi" = res_FRNA_mTORi[, "padj"],
                        
                        "baseMean_LRNA" = res_LRNA_2i[, "baseMean"],
                        "log2FoldChange_LRNA_2i" = res_LRNA_2i[, "log2FoldChange"],
                        "padj_LRNA_2i" = res_LRNA_2i[, "padj"],
                        "log2FoldChange_LRNA_mTORi" = res_LRNA_mTORi[, "log2FoldChange"],
                        "padj_LRNA_mTORi" = res_LRNA_mTORi[, "padj"],
                        
                        "baseMean_FRNA_sp" = res_FRNA_2i_sp[, "baseMean"],
                        "log2FoldChange_FRNA_sp_2i" = res_FRNA_2i_sp[, "log2FoldChange"],
                        "padj_FRNA_sp_2i" = res_FRNA_2i_sp[, "padj"],
                        "log2FoldChange_FRNA_sp_mTORi" = res_FRNA_mTORi_sp[, "log2FoldChange"],
                        "padj_FRNA_sp_mTORi" = res_FRNA_mTORi_sp[, "padj"],
                        
                        "baseMean_LRNA_sp" = res_LRNA_2i_sp[, "baseMean"],
                        "log2FoldChange_LRNA_sp_2i" = res_LRNA_2i_sp[, "log2FoldChange"],
                        "padj_LRNA_sp_2i" = res_LRNA_2i_sp[, "padj"],
                        "log2FoldChange_LRNA_sp_mTORi" = res_LRNA_mTORi_sp[, "log2FoldChange"],
                        "padj_LRNA_sp_mTORi" = res_LRNA_mTORi_sp[, "padj"]
  )
  
  if (.include_SL2i) {
    res_dat <- cbind(res_dat,
                     data.frame("log2FoldChange_FRNA_sp_2i" = res_FRNA_SL2i_sp[, "log2FoldChange"],
                                "padj_FRNA_sp_2i" = res_FRNA_SL2i_sp[, "padj"],
                                "log2FoldChange_LRNA_2i" = res_LRNA_SL2i[, "log2FoldChange"],
                                "padj_LRNA_2i" = res_LRNA_SL2i[, "padj"]))
  }
  
  mcols(TU.DE) <- res_dat
  TU.DE
}

TU.DE.mm10 <- get_DE_res(TU.counts.mat, TU.counts)
saveRDS(TU.DE.mm10, "data/TU.DE.mm10.RData")

TU.DE.mm9 <- get_DE_res(TU.counts.mat.mm9, TU.counts.mm9, .include_SL2i = T)
saveRDS(TU.DE.mm9, "data/TU.DE.mm9.RData")

# add attributes, direction, enhancer ----------------------------------------------------------------------------
if (T) {
  # enhancer in development
  F5_enhancer <- importRanges("../data/F5.mm10.enhancers.bed")
  
  TU.DE.intergenic <- TU.DE.mm10[TU.DE.mm10$location == "intergenic" & TU.DE.mm10$baseMean_LRNA > 0]
  TU.intergenic_as <- promoters(TU.DE.intergenic, upstream = 500, downstream = 0)
  levels(strand(TU.intergenic_as)) <- c("-", "+", "*")
  
  mtch_dir <- findOverlaps(TU.DE.intergenic, TU.intergenic_as)
  TU.DE.intergenic$direction <- ifelse(countQueryHits(mtch_dir) > 0,
                                       "Bidirectional", "Unidirectional" )
  TU.DE.intergenic$enhancer <- ifelse(findOverlaps(F5_enhancer, TU.DE.intergenic) %>% countSubjectHits > 0,
                                      "Enhancer", "TX")
  TU.DE.intergenic$enhancer_direction <- paste(TU.DE.intergenic$direction, TU.DE.intergenic$enhancer, sep = "_")
  # major / minor direction by expr
  mtch_dir <- as.matrix(mtch_dir)
  for( i in seq_len(nrow(mtch_dir)) )
  { # sort by reads density
    if (which.max(TU.DE.intergenic$baseMean_LRNA[mtch_dir[i, ]] /
                  width(TU.DE.intergenic[mtch_dir[i, ]])) == 2)
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

  phastConsnMat <- .countBW(phastCons.list[1], promoters(TU.DE.intergenic, 500, 200),
                            fast = F) # 60way is okay
  colnames(phastConsnMat) <- gsub(".*(mm10.*).bw", "\\1", colnames(phastConsnMat))
  mcols(TU.DE.intergenic) <- cbind(mcols(TU.DE.intergenic), phastConsnMat)
  
  TU.DE.intergenic$mm10.60way.phastCons <- phastConsnMat$mm10.60way.phastCons.bw
  
  # ChromHMM state coeffiencts on log2FC, uni- / bi-direction intergenic TU
  intergenic_2i_DE <- mcols(TU.DE.intergenic)[TU.DE.intergenic$enhancer == "TX", ] %>% as.data.frame()
  coef_dat <- rbind(
    cbind("Bidirection", summary(lm( log2FoldChange_LRNA_2i ~ chromHMM,
                                     intergenic_2i_DE[intergenic_2i_DE$direction == "Bidirectional", ]))$coef),
    cbind("Unidirection", summary(lm( log2FoldChange_LRNA_2i ~ chromHMM,
                                      intergenic_2i_DE[intergenic_2i_DE$direction == "Unidirectional", ]))$coef)
  ) %>% as.data.frame
  coef_dat$chromHMM <- gsub(".1", "\\2", gsub("(chromHMM|X.)", "\\2", rownames(coef_dat)))
  coef_dat$chromHMM[grep("Intercept", coef_dat$chromHMM)] <- "ActivePromoter"
  coef_dat$chromHMM <- factor(coef_dat$chromHMM, levels = sort(unique(coef_dat$chromHMM), T) )
  coef_dat <- coef_dat[!grepl("Transcription", coef_dat$chromHMM), ]
  coef_dat$Estimate <- as.numeric(as.character(coef_dat$Estimate))
  coef_dat$"Std. Error" <- as.numeric(as.character(coef_dat$"Std. Error"))
  saveRDS(coef_dat, "../fig2/data/coef_dat.RData")
  
  # Bidirectional intergenic TU
  bi_intergenic_DE <- with(TU.DE.intergenic[c(mtch_dir)],
                           data.frame(log2FoldChange = c(log2FoldChange_LRNA_2i, log2FoldChange_LRNA_mTORi),
                                      padj = c(padj_LRNA_2i, padj_LRNA_mTORi),
                                      Sample = c(rep("2i", length(mtch_dir)),
                                                 rep("mTORi", length(mtch_dir) ) ),
                                      Direction = c(rep("Major", nrow(mtch_dir)),
                                                    rep("Minor", nrow(mtch_dir))),
                                      Enhancer = rep(enhancer, 2) ))
  saveRDS(bi_intergenic_DE, "../fig2/data/bi_intergenic_DE.RData")
  
  # uni directional intergenic TU DE
  uni_intergenic_DE <- with(TU.DE.intergenic[-c(mtch_dir)],
                            data.frame(log2FoldChange = c(log2FoldChange_LRNA_2i, log2FoldChange_LRNA_mTORi),
                                       padj = c(padj_LRNA_2i, padj_LRNA_mTORi),
                                       Sample = c(rep("2i", length(TU.DE.intergenic) - length(mtch_dir)),
                                                  rep("mTORi", length(TU.DE.intergenic) - length(mtch_dir))),
                                       Enhancer = rep(enhancer, 2) ))
  saveRDS(uni_intergenic_DE, "../fig2/data/uni_intergenic_DE.RData")
  
  bi_intergenic.other.mm10.gr <- TU.DE.intergenic[intersect.Vector(c(mtch_dir),
                                                                   which(TU.DE.intergenic$enhancer == 'TX'))]
  bi_intergenic.other.mm9.gr <- liftOver(x = bi_intergenic.other.mm10.gr, 
                                         chain = import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>% Reduce(f = "c")
  
}

