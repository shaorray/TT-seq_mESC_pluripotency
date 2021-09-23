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
names(gene.gr) <- gene.gr$gene_id 

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
                                 stranded = T,
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

est_speed_mat <- TT_gene_density_mat[rownames(gene_body_mat), ] - as.matrix(gene_body_mat)
saveRDS(est_speed_mat, "../fig4/data/est_speed_mat.RData")

# combine Pol2 and RNA synthesis
gene_body_production_mat <- cbind(gene_body_mat, TT_gene_density_mat[rownames(gene_body_mat), ])

gene_body_production_mat$LRNA_RPK <- (TU.DE.mm9.gr$baseMean_LRNA_sp / width(TU.DE.mm9.gr) * 1e3)[
  match(rownames(gene_body_production_mat), TU.DE.mm9.gr$gene_id)] %>% log10()

gene_body_production_mat <- gene_body_production_mat[
  !apply(gene_body_production_mat, 1, function(x) any(is.infinite(as.numeric(x)) | is.na(x) | x == 0)), ]

# load 2014 Jonker et.al Fig3 data ------------------------------------------------------------------------
paper_est_speed <- read.table("../data/elife-02407-fig3-data1-v1", sep = "\t", header = 1)

gene_mtch <- biomaRt::getBM(filters = "refseq_mrna", 
                            attributes = c("refseq_mrna", "ensembl_gene_id"),
                            values = paper_est_speed$Genename, 
                            mart = biomaRt::useMart(dataset = "mmusculus_gene_ensembl", 
                                                    biomart = "ensembl"))

paper_est_speed <- paper_est_speed[match(gene_mtch$refseq_mrna, paper_est_speed$Genename), ]
paper_est_speed$gene_id <- gene_mtch$ensembl_gene_id
paper_est_speed$speed <- paper_est_speed$Rate..bp.min. / 1e3

######################################################################################################
# Pol II pausing index, start-seq TSS interval / (+500, +1500) gene body
######################################################################################################
Start_peaks.gr <- reduce(importRanges("../data/Start_seq_peaks.gtf"), min.gapwidth=50L)
tss.gr <- promoters(gene.gr, upstream = 0, downstream = 0)

gene_bodies.gr <- flank(promoters(tss.gr, upstream = 0, downstream = 500), width = 1500, start = F)

start_hits <- findOverlaps(Start_peaks.gr, tss.gr+100, ignore.strand = F) # adjust TSS to the nearest Start-seq peak
pause_sites.gr <- GRanges()
for(i in seq_along(tss.gr))
{
  tmp_hits <- queryHits(start_hits)[subjectHits(start_hits) == i]
  if (length(tmp_hits) == 1) {
    pause_sites.gr <- c(pause_sites.gr, Start_peaks.gr[tmp_hits])
  } else if (length(tmp_hits) > 1) {
    keep.idx <- nearest(tss.gr[i], Start_peaks.gr[tmp_hits], ignore.strand = F)
    pause_sites.gr <- c(pause_sites.gr, Start_peaks.gr[tmp_hits][keep.idx] )
  } else {
    pause_sites.gr <- c(pause_sites.gr, promoters(tss.gr[i], upstream = 0, downstream = 50))
  }
}
names(pause_sites.gr) <- NULL

if (T) {
  # MINUTE ChIP data
  pausing_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1.*Pol5p.*ALL.*bam$", full.names = T),
                           intervals = pause_sites.gr, stranded = F, paired.end = "ignore")
  pausing_mat <- pausing_mat / width(pause_sites.gr) # convert to reads density
  
  gene_body_mat <- .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1.*Pol5p.*ALL.*bam$", full.names = T),
                             intervals = gene_bodies.gr, stranded = F, paired.end = "ignore")
  gene_body_mat <- gene_body_mat / width(gene_bodies.gr) # convert to reads density
  
  pausing_index_mat <- pausing_mat / gene_body_mat
}
# pausing_index_mat <- pausing_index_mat[!is.na(tss.gr$gene_id), ]
rownames(pausing_index_mat) <- tss.gr$gene_id
saveRDS(pausing_index_mat, "data/pausing_index_SL_2i.RData")


# ------------------------------------- Review supports --------------------------------------------- #
# compare est velocity by different Pol2 ChIP 
if (F) {
  bw_files <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/bw_scaled",
                         pattern = "Pol.*(S5|S2).*pool", full.names = T)
  bw_files <- bw_files[!grepl("2iSL|2i_INH", bw_files)]
  
  Pol_bw_gene_body_mat <- .countBW(bw_files = bw_files, intervals = gene.gr, fast = F)
  Pol_bw_gene_body_mat <- Pol_bw_gene_body_mat[, c(2,1,3, 5,4,6)] %>% log10()
  
  Pol_est_speed_cmp <- cbind(TT_gene_density_mat - Pol_bw_gene_body_mat[, c(2,1,3)] %>% "/"(width(gene.gr)) %>% log10(),
                             TT_gene_density_mat - Pol_bw_gene_body_mat[, c(5,4,6)] %>% "/"(width(gene.gr)) %>% log10()) %>% 
    trim_quantile() %>% as.data.frame()
  colnames(Pol_est_speed_cmp) <- c("SL_Pol2S2p", "2i_Pol2S2p", "mTORi_Pol2S2p", "SL_Pol2S5p", "2i_Pol2S5p", "mTORi_Pol2S5p")
  
  plot_scatter <- function(dat, .xlab, .ylab, xlim = NULL, ylim = NULL) {
    dat <- dat[complete.cases(dat) & is.finite(rowSums(dat)) & !is.na(rowSums(dat)), ]
    r <- single_variance_explained(dat[, 1], dat[, 2], T) %>% round(3)
    ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
      geom_point(cex = 0.5) +
      annotate("text", x = -Inf, y = Inf,
               hjust = -0.5, vjust = 1.2, 
               label = paste0(" r = ", r, "\nn = ", nrow(dat))) +
      scale_x_continuous(name = .xlab, limits = xlim) +
      scale_y_continuous(name = .ylab, limits = ylim) +
      scale_color_viridis_c(option = "A", direction = -1, begin = 0.1, end = 0.9) +
      theme_setting +
      theme(legend.position = "none")
  }
  
  g1 <- data.frame(x = Pol_est_speed_cmp$SL_Pol2S2p, 
                   y = Pol_est_speed_cmp$SL_Pol2S5p) %>% 
    plot_scatter(.xlab = "Est. velocity (TTseq / PolIIS2p) (a.u.)", 
                 .ylab = "Est. velocity (TTseq / PolIIS5p) (a.u.)",
                 r = "0.991") + ggtitle("SL")
  
  g2 <- data.frame(x = Pol_est_speed_cmp$`2i_Pol2S2p`, 
                   y = Pol_est_speed_cmp$`2i_Pol2S5p`) %>% 
    plot_scatter(.xlab = "Est. velocity (TTseq / PolIIS2p) (a.u.)", 
                 .ylab = "Est. velocity (TTseq / PolIIS5p) (a.u.)",
                 r = "0.988") + ggtitle("2i")
  
  g3 <- data.frame(x = Pol_est_speed_cmp$mTORi_Pol2S2p, 
                   y = Pol_est_speed_cmp$mTORi_Pol2S5p) %>% 
    plot_scatter(.xlab = "Est. velocity (TTseq / PolIIS2p) (a.u.)", 
                 .ylab = "Est. velocity (TTseq / PolIIS5p) (a.u.)",
                 r = "0.979") + ggtitle("mTORi")
  
  ggsave(grid.arrange(g1, g2, g3, nrow = 1), 
         filename = "../figS3/figs/FigS3_scatter_est_velo_Pol2_types.png",
         width = 11.1, height = 3.7)
  
  # NET-seq
  bam_files = c("/mnt/0E471D453D8EE463/GEO_nascent_RNA_mm9/2018_Mylonas_NET-seq_Control_rep1/bam/SRR5350576_pass.mm9.Aligned.sortedByCoord.out.bam", 
                "/mnt/0E471D453D8EE463/GEO_nascent_RNA_mm9/2018_Mylonas_NET-seq_Control_rep2/bam/SRR5350577_pass.mm9.Aligned.sortedByCoord.out.bam")
  
  gencode.mm9 <-
    importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
  gencode.mm9$gene_id <- gsub("\\..*", "", gencode.mm9$gene_id)
  gencode.mm9 <- gencode.mm9[gencode.mm9$gene_id %in% names(gene.gr) & gencode.mm9$type == "exon"]
  
  
  NET_exon_count <- Rsubread::featureCounts(
    bam_files,
    annot.ext = data.frame(
      GeneID = gencode.mm9$gene_id,
      Chr = as.character(seqnames(gencode.mm9)),
      Start = start(gencode.mm9),
      End = end(gencode.mm9),
      Strand = strand(gencode.mm9)
    ),
    isPairedEnd = F,
    strandSpecific = 1,
    nthreads = 4
  )
  NET_exon_RPK <- rowMeans(NET_exon_count$counts) / NET_exon_count$annotation$Length * 1e3
  
  TT_exon_count <- Rsubread::featureCounts(
    c("/mnt/0E471D453D8EE463/TT_seq_data/RS20190401-125485360/bam_mm9/LRNA_SL_rep1.Aligned.sortedByCoord.out.bam", 
      "/mnt/0E471D453D8EE463/TT_seq_data/RS20190401-125485360/bam_mm9/LRNA_SL_rep2.Aligned.sortedByCoord.out.bam"),
    annot.ext = data.frame(
      GeneID = gencode.mm9$gene_id,
      Chr = as.character(seqnames(gencode.mm9)),
      Start = start(gencode.mm9),
      End = end(gencode.mm9),
      Strand = strand(gencode.mm9)
    ),
    isPairedEnd = T,
    strandSpecific = 1,
    nthreads = 4
  )
  TT_exon_RPK <- rowMeans(TT_exon_count$counts) / TT_exon_count$annotation$Length * 1e3
  
  g4 <- data.frame(x = TT_gene_density_mat[rownames(est_speed_mat), 1] - log10(NET_exon_RPK[rownames(est_speed_mat)]), 
                   y = est_speed_mat[, 1]) %>% 
    plot_scatter(.xlab = "Est. velocity (TTseq / NET-seq) (a.u.)", 
                 .ylab = "Est. velocity (TTseq / PolS5p) (a.u.)") 
  
  # pausing estimation
  NET_pause_count <- .countBam(bam_files, intervals = pause_sites.gr,
                               stranded = T, paired.end = "ignore")
  
  TT_puase_count <- .countBam(bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                     pattern = "LRNA.*(SL_).*bam$", full.names = T), 
                              intervals = pause_sites.gr,
                              stranded = T, paired.end = "ignore")
  
  Pol2S5p_pause_count <- .countBam(bam_files = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/P1_Pol5p_SL_CTR_ALL.mm9.fltd.bam", 
                                   intervals = pause_sites.gr,
                                   stranded = F, paired.end = "ignore")
  
  NET_gene_count <- .countBam(bam_files, intervals = gene.gr,
                               stranded = T, paired.end = "ignore")
  Pol2S5p_gene_count <- .countBam(bam_files = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/P1_Pol5p_SL_CTR_ALL.mm9.fltd.bam", 
                                   intervals = gene.gr,
                                   stranded = F, paired.end = "ignore")
  
  
  puase_dat <- data.frame(x = pausing.time[, 1],
                          y = log10(width(pause_sites.gr) / 1e3) -
                                log10(rowMeans(TT_puase_count) + 1) +
                                log10(rowMeans(NET_pause_count) + 1))
  rm.idx <- which(pausing.time[, 1] %in% names(head(sort(table(pausing.time[,1]), T))))
  
  g5 <- plot_scatter(puase_dat[-rm.idx, ], .xlab = "Est.pausing time(TTseq/NETseq)(a.u.)", 
                 .ylab = "Est. pausing time (TTseq/PolS5p) (a.u.)")
  
  # NET-seq speed vs measured speed
  gene_intersect <- intersect.Vector(rownames(elongation_speed_table), rownames(TT_gene_density_mat))
  gene_intersect <- gene_intersect[elongation_speed_table[gene_intersect, 1] > 0.1]
  g6 <- data.frame(x = (TT_gene_density_mat[gene_intersect, 1] - log10(NET_exon_RPK[gene_intersect])) %>% trim_quantile(), 
             y = elongation_speed_table[gene_intersect, 1]) %>% 
    plot_scatter(.xlab = "Est. velocity (TTseq / NET-seq) (a.u.)", 
                 .ylab = "") + 
    scale_y_log10(name = "Measured velocity (kb/min)")

  ggsave(grid.arrange(g4, g5, g6, nrow = 1), 
         filename = "../figS3/figs/FigS3_scatter_est_pausing_time_NET.png",
         width = 12, height = 4)
}

if (F) {
  # correlation between ChIP replicates
  bam_files <- list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam", "P1_Pol5p.*_R.*bam$", full.names = T)
  Pol5p_reps <- .countBam(bam_files = bam_files, intervals = gene.gr, stranded = F, paired.end = "ignore")
  Pol5p_reps <- log(Pol5p_reps / width(gene.gr) * 1e3) # log RPK without normalization
  
  g1.1 <- plot_scatter(`colnames<-`(Pol5p_reps[, 7:8], c("x", "y")), 
                       .xlab = "SL rep1", .ylab = "SL rep2")
  g1.2 <- plot_scatter(`colnames<-`(Pol5p_reps[, 8:9], c("x", "y")), 
                       .xlab = "SL rep2", .ylab = "SL rep3")
  
  g2.1 <- plot_scatter(`colnames<-`(Pol5p_reps[, 1:2], c("x", "y")), 
               .xlab = "2i rep1", .ylab = "2i rep2")
  g2.2 <- plot_scatter(`colnames<-`(Pol5p_reps[, 2:3], c("x", "y")), 
                       .xlab = "2i rep2", .ylab = "2i rep3")
  
  g3.1 <- plot_scatter(`colnames<-`(Pol5p_reps[, 10:11], c("x", "y")), 
                       .xlab = "mTORi rep1", .ylab = "mTORi rep2")
  g3.2 <- plot_scatter(`colnames<-`(Pol5p_reps[, 11:12], c("x", "y")), 
                       .xlab = "mTORi rep2", .ylab = "mTORi rep3")

  ggsave(grid.arrange(g1.1, g2.1, g3.1, g1.2, g2.2, g3.2, nrow = 2), 
         filename = "../figS3/figs/FigS3_scatter_Pol2S5p_replicates.png",
         width = 9, height = 4)
}

# ---------------------------- SL2i batch ---------------------------------- #
# read Pol2s5p gene body + flanks (sandwich) coverage
bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                       pattern = "P2.*Pol.*ALL.*bam$", full.names = T)

pol2s5p_sl2i_mat_sandwich <- readBam(bam_files = bam_files,
                                intervals = gene.gr,
                                pair_end = T,
                                stranded = F,
                                flanks = c(2000, 4000),
                                new_lens = c(20, 200, 40))

pol2s5p_sl2i_input_size_factor = .countBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                                             pattern = "P2_IN_.*ALL.*bam$", full.names = T)[-1],
                                      gene.gr) %>% SizeFactorCal()

pol2s5p_sl2i_mat_sandwich_list <- list()
for (i in c(1, 3)) {
  pol2s5p_sl2i_mat_sandwich_list <- c(pol2s5p_sl2i_mat_sandwich,
                                 list(pol2s5p_sl2i_mat_sandwich[[i]] / pol2s5p_sl2i_input_size_factor[i]))
}
names(pol2s5p_sl2i_mat_sandwich_list) <- c("SL", "SL2i_2d")

# TTseq LRNA coverage
bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                       pattern = "LRNA.*(SL_|_SL2i_2d).*bam$", full.names = T)

TT_sl2i_gene_sandwich_mat <- readBam(bam_files = bam_files,
                                     intervals = gene.gr,
                                     pair_end = T,
                                     stranded = F, # include upstream reads on the main gene coverage
                                     flanks = c(2000, 4000),
                                     new_lens = c(20, 200, 40))

names(TT_sl2i_gene_sandwich_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1", names(TT_sl2i_gene_sandwich_mat))

# normalise TTseq coverage by spike-in size factors
for (i in names(TT_sl2i_gene_sandwich_mat)) {
  TT_sl2i_gene_sandwich_mat[[i]] <- TT_sl2i_gene_sandwich_mat[[i]] / L_RNA_sf[i]
}
# averaging replicates
TT_sl2i_gene_sandwich_mat_norm <- list()
for (i in c("SL_", "SL2i_2d")) {
  idx <- grep(i, names(TT_sl2i_gene_sandwich_mat))
  TT_sl2i_gene_sandwich_mat_norm <- c(TT_sl2i_gene_sandwich_mat_norm,
                                 list(Reduce("+", TT_sl2i_gene_sandwich_mat[idx]) / length(idx)))
}

names(TT_sl2i_gene_sandwich_mat_norm) <- c("SL", "SL2i_2d")

# ---
TT_Pol2s5p_sl2i_sandwich_mat <- list()
for (i in c("SL", "SL2i_2d")) {
  TT_Pol2s5p_sl2i_sandwich_mat <- c(TT_Pol2s5p_sl2i_sandwich_mat,
                               list(log1p(TT_sl2i_gene_sandwich_mat_norm[[i]]) - log1p(pol2s5p_sl2i_mat_sandwich_list[[i]]) * 5) )
}
names(TT_Pol2s5p_sl2i_sandwich_mat) <- c("SL", "SL2i_2d")


if (F) {
  # Pol2S5p
  Pol2S5p_SL2i_gene_body_mat <- 
    .countBam(bam_files = list.files("/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                     "P2.*Pol5p.*ALL.*bam$", full.names = T)[c(1, 3)],
              intervals = gene.gr,
              stranded = F, paired.end = "ignore")
  
  Pol2S5p_SL2i_gene_body_mat <- sweep(Pol2S5p_SL2i_gene_body_mat, 2,
                                      pol2s5p_sl2i_input_size_factor[c(1, 3)], "/")
  # TT-seq
  bam_files = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                         pattern = "LRNA.*(SL_|SL2i_2d).*(rep1|rep2).*bam$", full.names = T)
  
  TT_SL2i_gene_density_mat <- .countBam(bam_files = bam_files,
                                        intervals = gene.gr,
                                        stranded = T,
                                        paired.end = "ignore")
  
  names(TT_SL2i_gene_density_mat) <- gsub("LRNA_(.*).Aligned.*", "\\1",
                                          names(TT_SL2i_gene_density_mat))
  
  # normalise TT-seq coverage by spike-in size factors
  TT_SL2i_gene_density_mat <- sweep(TT_SL2i_gene_density_mat, 2, 
                                    L_RNA_sf[names(TT_SL2i_gene_density_mat)], "/") 
  
  # averaging replicates
  tmp <- NULL
  for (i in c("SL_", "SL2i_")) {
    idx <- grep(i, colnames(TT_SL2i_gene_density_mat))
    tmp <- cbind(tmp,
                 TT_SL2i_gene_density_mat[, idx] %>%
                   "/"(width(gene.gr) / 1e3) %>% 
                   rowMeans() %>% log10() )
  }
  TT_SL2i_gene_density_mat <- `rownames<-`(tmp, gene.gr$gene_id)
  colnames(TT_SL2i_gene_density_mat) <- c("Tx_RPK_SL", "Tx_RPK_SL2i")
  
  est_SL2i_speed_mat <- 
    TT_SL2i_gene_density_mat - 
    as.matrix(log10(Pol2S5p_SL2i_gene_body_mat))
  est_SL2i_speed_mat <- est_SL2i_speed_mat[!is.na(rowSums(est_SL2i_speed_mat)), ]
  
  # "mRNA.speed" is from "Fig4_Est_speed_explain.R" line 112
  # scaling est. velocity
  genes_ov <- intersect.Vector(rownames(est_SL2i_speed_mat), rownames(mRNA.speed))
  est_SL2i_speed_mat <-
    est_SL2i_speed_mat[genes_ov,] - 
    median(est_SL2i_speed_mat[genes_ov, 1] - mRNA.speed[genes_ov, 1])
  
  est_SL2i_speed_mat <- cbind(mRNA.speed[genes_ov, ], "SL2i" = est_SL2i_speed_mat[, 2])
  est_SL2i_speed_mat <- est_SL2i_speed_mat[is.finite(rowSums(est_SL2i_speed_mat)), ]
  
  # plot sample est. velocity
  data_speed_all <- reshape::melt(est_SL2i_speed_mat)
  data_speed_all$X2 <- factor(data_speed_all$X2, c("SL", "2i", "SL2i", "mTORi"))
  
  ggplot(data_speed_all, aes(x = X2, y = value, fill = X2)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(name = "Est. velocity (kb/min)", 
                       breaks = c((-2):2), labels = 10^c((-2):2)) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    ggpubr::stat_compare_means(label = "p.signif", method = "t.test", 
                               ref.group = "SL") +
    scale_fill_manual(values = colors_20[c(13, 2, 10, 7)]) + 
    xlab("") + ggtitle("mRNA (n = 10617)") +
    theme_setting +
    theme(legend.position = "none")
  
  ggsave(filename = "FigS3_Est_speed_mRNA_boxplot.png",
         path = "../figS3/figs",
         device = "png", width = 4, height = 4)
  
}
