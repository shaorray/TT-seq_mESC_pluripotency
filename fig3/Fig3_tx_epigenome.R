#-----------------------------------------------------------------------------------------------------
# This part links epigenome states, regulatory annotations to TU annotation's tx activity, half-life, 
# copy per cell, in the aim to understand the internal changes took places in pluripotent states
#
# Rui Shao, Jul 2020
#-----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# Intergenic TU ChromHMM states --------------------------------------------------------------------------------------------
sample_ncRNA_counts_Rates <- readRDS("../figS2/data/sample_nc_counts_Rates_combined.RData")
ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub("(.*?)_(.*)", "\\2", ChromHMM$name)
# ChromHMM <- liftOver(ChromHMM,
#                      chain = import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")
# ) %>% unlist()

TU.gr <- sample_ncRNA_counts_Rates %>%
  dplyr::filter(!grepl("sense", gene_id) & Sample == "SL") %>%
  "$"(gene_id) %>%
  as.character() %>%
  strsplit("_") %>%
  {Reduce(rbind, .)} %>%
  as.data.frame() %>%
  `colnames<-`(c("location", "seqname", "start", "end", "strand")) %>%
  makeGRangesFromDataFrame()

TU.gr$location <- with(sample_ncRNA_counts_Rates,
                       gsub("_.*", "\\1", gene_id[!grepl("sense", gene_id) & Sample == "SL"]))
# TSS state
TU.gr$chromHMM <- NA
mtch <- findOverlaps(promoters(TU.gr, 1, 0), ChromHMM)
TU.gr$chromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]
# Gene body ChromHMM longest state
TU.gr$chromHMM_gb <- NA
mtch <- findOverlaps(TU.gr, ChromHMM, minoverlap = 200)
for (i in unique(queryHits(mtch))) {
  tmp.ChromHMM <- ChromHMM[subjectHits(mtch)[queryHits(mtch) == i]]
  tot_w <- aggregate(width(tmp.ChromHMM), list(tmp.ChromHMM$name), "sum")
  TU.gr$chromHMM_gb[i] <- tot_w[which.max(tot_w[, "x"]), 1]
}
# Direction
TU.as <- promoters(TU.gr, upstream = 1000, downstream = 500)
levels(strand(TU.as)) <- c("-", "+", "*")
mtch_dir <- findOverlaps(TU.gr, TU.as)
TU.gr$Direction <- ifelse(mtch_dir %>% countQueryHits() > 0,
                          "Bidirection", "Unidirection" )
# CpG islands
TU.gr$CpG_counts <- interval_pattern_hits(intervals = promoters(TU.gr, 
                                                                upstream = 1000,
                                                                downstream = 50),
                                          pattern = "CG", 
                                          which_genome = "mm10")
TU.gr$CpG_islands <- TU.gr$CpG_counts > sort(TU.gr$CpG_counts, T)[which.max(cumsum( sort(scale(TU.gr$CpG_counts), T) ))]

# Heatmap ChromHMM groups ----------------------------------------------------------------
dat6 <- sample_ncRNA_counts_Rates %>%
  {as.data.frame(split(.$Copy, .$Sample))} %>%
  `colnames<-`(unique(sample_ncRNA_counts_Rates$Sample)) %>%
  dplyr::select(c("SL", "2i_2d", "mTORi_1d")) %>%
  `rownames<-`(unique(sample_ncRNA_counts_Rates$gene_id)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(contains("intergenic")) %>%
  t()

dat7 <- sample_ncRNA_counts_Rates %>%
  {as.data.frame(split(.$Half_life, .$Sample))} %>%
  `colnames<-`(unique(sample_ncRNA_counts_Rates$Sample)) %>%
  dplyr::select(c("SL", "2i_2d", "mTORi_1d")) %>%
  `rownames<-`(unique(sample_ncRNA_counts_Rates$gene_id)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(contains("intergenic")) %>%
  t()

annotation_row <- data.frame(TSS.state = TU.gr$chromHMM[TU.gr$location == "intergenic"],
                             Body.state = TU.gr$chromHMM_gb[TU.gr$location == "intergenic"],
                             Direction = TU.gr$Direction[TU.gr$location == "intergenic"],
                             CpG.island = c("Low", "High")[TU.gr$CpG_islands[TU.gr$location == "intergenic"] + 1],
                             row.names = rownames(dat7))

dat8 <- cbind((log10(dat6)), (log10(dat7)))[order(annotation_row$TSS.state), ]
dat8 <- dat8[!apply(dat8, 1, function(x) any(is.na(x) | is.infinite(x) | x < -2)), ]

set.seed(1)
cls <- kmeans(dat8, centers = 4)$cluster
new_dat8 <- dat8[order(cls), ]
row_state <- annotation_row[rownames(new_dat8), "TSS.state"]
cls <- sort(cls)

out_dat <- NULL
for (i in c(1,2,3,4)) {
  # each cluster
  tmp_state <- row_state[cls == i]
  tmp_dat <- new_dat8[cls == i, ][order(tmp_state), ]
  tmp_state <- sort(tmp_state)
  for (st in unique(tmp_state)) {
    # each state
    if (sum(tmp_state == st) > 1) {
      out_dat <- rbind(out_dat, 
                       tmp_dat[tmp_state == st, ][order(rowSums(tmp_dat[tmp_state == st, c(1:3)]), decreasing = T), ])
    } else {
      out_dat <- rbind(out_dat, tmp_dat[tmp_state == st, ])
    }
  }
}

out_dat[out_dat > 3] <- 3; out_dat[out_dat < (-1)] <- -1; out_dat <- out_dat[nchar(rownames(out_dat)) > 0, ]
annotation_row <- annotation_row[rownames(out_dat), ]
annotation_colors <- list(TSS.state = colors_20[c(1,2,3,5,7,8,9,11,12,13,14)],
                          Body.state = colors_20[c(1,2,3,5,7,8,9,11,12,13,14)],
                          CpG.island = colors_9[4:3],
                          Direction = colors_9[1:2])
names(annotation_colors[[1]]) <- unique(row_state)[-12]
names(annotation_colors[[2]]) <- unique(row_state)[-12]
names(annotation_colors[[3]]) <- c("High", "Low")
names(annotation_colors[[4]]) <- c("Bidirection", "Unidirection")
pdf("figs/Fig3_Intergenic_copy_half_life_ChromHMM.pdf", width = 5, height = 4.5)
out <- pheatmap::pheatmap(out_dat,
                          color = rev(viridis::inferno(100, begin = 0.2)),
                          border_color = NA, cellwidth = NA, cellheight = NA,
                          legend_breaks = (-1):3,
                          legend_labels = c("<0.1", 1, 10, 100, ">1000"),
                          angle_col = 315,
                          labels_col = c("SL", "2i 2d", "mTORi 1d",
                                         "SL", "2i 2d", "mTORi 1d"),
                          cluster_cols = F,
                          cluster_rows = F,
                          show_rownames = F,
                          annotation_row = annotation_row,
                          annotation_colors = annotation_colors,
                          gaps_row = cumsum(table(cls[rownames(out_dat)])[c(1,2,3,4)]),
                          gaps_col = 3)
dev.off()

# TX ChromHMM states --------------------------------------------------------------------------------------------
# ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
# ChromHMM$name <- gsub("(.*?)_(.*)", "\\2", ChromHMM$name)
ChromHMM <- liftOver(ChromHMM,
                     chain = import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")
) %>% unlist()

tss.gr <- importRanges("../data/tss.attributes.gff3")
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(tss.gr$gene_name),
                       keytype = "GENENAME",
                       columns = c("ENTREZID", "GENEID"))
tss.gr$gene_id <- res$GENEID[match(tss.gr$gene_name, res$GENENAME)]
tx.tss.gr <- tss.gr[!is.na(tss.gr$gene_id)]

sample_tx_counts_Rates <- readRDS("../figS2/data/sample_Tx_counts_Rates_combined.RData")

dat3 <- sample_tx_counts_Rates %>%
  {as.data.frame(split(.$Copy, .$Sample))} %>%
  `colnames<-`(unique(sample_tx_counts_Rates$Sample)) %>%
  dplyr::select(c("SL", "2i_2d", "mTORi_1d")) %>%
  t() %>%
  as.data.frame() %>%
  t()

dat4 <- sample_tx_counts_Rates %>%
  {as.data.frame(split(.$Half_life, .$Sample))} %>%
  `colnames<-`(unique(sample_tx_counts_Rates$Sample)) %>%
  dplyr::select(c("SL", "2i_2d", "mTORi_1d")) %>%
  t() %>%
  as.data.frame() %>%
  t()

dat5 <- cbind(dat3, dat4) %>% log10
rownames(dat5) <- unique(sample_tx_counts_Rates$gene_id)
dat5 <- dat5[match(tx.tss.gr$gene_id, rownames(dat5)), ]
dat5 <- dat5[apply(dat5, 1, function(x) all(!is.na(x) & is.finite(x) & x > 0)), ]

tx.tss.gr <- tx.tss.gr[match(rownames(dat5), tx.tss.gr$gene_id)]
# TSS ChromHMM state
tx.tss.gr$chromHMM <- NA
mtch <- findOverlaps(promoters(tx.tss.gr, 1, 0), ChromHMM)
tx.tss.gr$chromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]

# Gene body ChromHMM longest state
gene.gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)[tx.tss.gr$gene_id]
gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub("(.*?)_(.*)", "\\2", ChromHMM$name)

tx.tss.gr$chromHMM_gb <- NA
mtch <- findOverlaps(gene.gr, ChromHMM, minoverlap = 200)
for (i in unique(queryHits(mtch))) {
  tmp.ChromHMM <- ChromHMM[subjectHits(mtch)[queryHits(mtch) == i]]
  tot_w <- aggregate(width(tmp.ChromHMM), list(tmp.ChromHMM$name), "sum")
  tx.tss.gr$chromHMM_gb[i] <- tot_w[which.max(tot_w[, "x"]), 1]
}

# Directionality
gene.gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
gene.gr <- gene.gr[gene.gr$gene_biotype == "protein_coding"]
promoter.gr <- promoters(gene.gr, upstream = 2000, downstream = 0)
levels(strand(promoter.gr)) <- c("-", "+", "*")
tx.tss.gr$Direction <- "Unidirection"
bi_idx <- findOverlaps(tx.tss.gr, c(TU.as[TU.as$location == "uaRNA"], promoter.gr)) %>% countQueryHits() > 0
tx.tss.gr$Direction[bi_idx] <- "Bidirection"
names(tx.tss.gr) <- tx.tss.gr$gene_id
rm(bi_idx)

# CpG islands
tx.tss.gr$CpG_counts <- interval_pattern_hits(intervals = promoters(tx.tss.gr, 
                                                                upstream = 1000,
                                                                downstream = 50),
                                          pattern = "CG", 
                                          which_genome = "mm9")
tx.tss.gr$CpG_islands <- tx.tss.gr$CpG_counts > sort(tx.tss.gr$CpG_counts, T)[which.max(cumsum( sort(scale(tx.tss.gr$CpG_counts), T) ))]

# clusterin on half-life and copy number
cls <- kmeans(dat5, centers = 4)$cluster
# order by clusters
new_dat5 <- dat5[order(cls), ]
row_state <- setNames(tx.tss.gr$chromHMM, tx.tss.gr$gene_id)[order(cls)]
cls <- sort(cls)

out_dat <- NULL
for (i in c(4,2,1,3)) {
  # each cluster
  tmp_state <- row_state[cls == i]
  tmp_dat <- new_dat5[cls == i, ][order(tmp_state), ]
  tmp_state <- sort(tmp_state)
  for (st in unique(tmp_state)) {
    # each state
    if (sum(tmp_state == st) > 1) {
      out_dat <- rbind(out_dat, 
                       tmp_dat[tmp_state == st, ][order(rowSums(tmp_dat[tmp_state == st, c(1:3)]), decreasing = T), ])
    } else {
      out_dat <- rbind(out_dat, tmp_dat[tmp_state == st, ])
    }
  }
}

out_dat[out_dat > 3] <- 3; out_dat[out_dat < (-1)] <- -1; out_dat <- out_dat[nchar(rownames(out_dat)) > 0, ]
annotation_row <- data.frame(TSS.state = row_state[rownames(out_dat)],
                             Body.state = tx.tss.gr[rownames(out_dat)]$chromHMM_gb,
                             Direction = tx.tss.gr[rownames(out_dat)]$Direction,
                             CpG.island = c("Low", "High")[tx.tss.gr[rownames(out_dat)]$CpG_islands + 1],
                             row.names = rownames(out_dat))

annotation_colors <- list(Body.state = colors_20[c(1:3,5,7:9,11:14)],
                          TSS.state = colors_20[c(1:3,5,7:9,11:14)],
                          CpG.island = colors_9[4:3],
                          Direction = colors_9[1:2])
names(annotation_colors[[1]]) <- unique(row_state)[-12]
names(annotation_colors[[2]]) <- unique(row_state)[-12]
names(annotation_colors[[3]]) <- c("High", "Low")
names(annotation_colors[[4]]) <- c("Bidirection", "Unidirection")

pdf("figs/Fig3_mRNA_copy_half_life_ChromHMM.pdf", width = 5, height = 4.5)
out <- pheatmap::pheatmap(out_dat,
                          color = rev(viridis::inferno(100, begin = 0.2)),
                          border_color = NA, cellwidth = NA, cellheight = NA,
                          legend_breaks = (-1):3,
                          legend_labels = c("<0.1", 1, 10, 100, ">1000"),
                          angle_col = 315,
                          labels_col = c("SL", "2i 2d", "mTORi 1d",
                                         "SL", "2i 2d", "mTORi 1d"),
                          cluster_cols = F,
                          cluster_rows = F,
                          show_rownames = F,
                          annotation_row = annotation_row,
                          annotation_colors = annotation_colors,
                          gaps_row = cumsum(table(cls[rownames(out_dat)])[c(4,2,1,3)]),
                          gaps_col = 3)
dev.off()

# non coding TU: intergenic differential expression ----------------------------------------
readRDS("data/bi_intergenic_DE.RData") %>% # from "F1_src_TU_anno.R"
  ggplot(aes(x = interaction(Enhancer, Sample),
             y = log2FoldChange, fill = Sample)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_boxplot(aes(linetype = Direction), outlier.size = 0, notch=TRUE) +
  ylim(c(-9, 9)) +
  scale_linetype_manual(name = "Direction", values = c(1, 11)) +
  scale_x_discrete(labels = c("2i \nEnhancer", "2i\nOther", "mTORi\nEnhancer", "mTORi\nOther")) +
  scale_fill_manual(values = colors_20[c(2, 20)]) +
  annotate("text", x=1.4, y=9,label = "p_value < 1e-5", cex=4,parse=TRUE) +
  annotate("text", x= c(t(matrix(c(1:4 - .18, 1:4 + .18), nrow = 4))), y=8,
           label= c("ns", "ns", "'*'", "'*'", "ns", "ns", "ns", "ns"), cex=4, parse=TRUE) +
  xlab("") +
  labs(title="Bidirectional Intergenic TU",x="") +
  theme_setting -> g3

readRDS("data/uni_intergenic_DE.RData") %>%
  ggplot(aes(x = interaction(Enhancer, Sample),
             y = log2FoldChange, fill = Sample)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_boxplot(outlier.size = 0, width = 0.5, notch=TRUE) +
  ylim(c(-9, 9)) +
  scale_fill_manual(values = colors_20[c(2, 20)]) +
  scale_x_discrete(labels = c("Enhancer", "Other", "mTORi\nEnhancer", "mTORi\nOther")) +
  annotate("text", x=1.5, y= 9,label = "p_value < 1e-5", cex=4,parse=TRUE) +
  annotate("text", x= 1:4, y= 8, label= rep("ns", 4), cex=4, parse=TRUE) +
  xlab("") +
  labs(title="Unidirectional Intergenic TU",x="") +
  theme_setting +
  theme(legend.position = "none") -> g4

# ChromHMM states overlap coefficient
coef_dat <- readRDS("data/coef_dat.RData")
ggplot(coef_dat, aes(x = chromHMM, y = Estimate, color = V1)) +
  geom_hline(yintercept = 0, linetype="solid", color = "grey", size= .5) +
  geom_errorbar(aes(ymin = Estimate - coef_dat$`Std. Error`,
                    ymax = Estimate + coef_dat$`Std. Error`), width = .2,
                position = position_dodge(.5)) +
  geom_point(position = position_dodge(0.5)) +
  scale_color_manual(values = colors_9[c(6, 2)]) +
  coord_flip() +
  labs(color="", title = "ChromHMM") +
  xlab("") + ylab("log2FC coef.") +
  theme_setting -> g5

ggsave(plot = egg::ggarrange(g4, g3, g5,nrow = 1, widths = c(6,6,4)),
       filename = "Fig3_non_coding_TU_DE.pdf", path = "figs",
       device = "pdf", width = 12, height = 4 )

# asRNA and intergenic TU overlapping with ChromHMM states --------------------------------
if (F)
{
  ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
  
  pacman::p_load(ggalluvial)
  
  nc_TU <- TU_SL_LRNA[TU_SL_LRNA$type == "ncRNA"]
  nc_TU <- nc_TU[!nc_TU$location %in% c("sense_intragenic", "dsRNA") ]
  nc_TU$location[nc_TU$location == "antisense"] <- "asRNA"
  nc_TU$ChromHMM <- NA
  
  chrm_TUs <- data.frame(Type = rep(unique(nc_TU$location), each = length(unique(ChromHMM$name))),
                         State = rep(unique(ChromHMM$name), length(unique(nc_TU$location))),
                         Hits = 0)
  
  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  
  for( type in unique(nc_TU$location) )
  {
    tmp_TU <- nc_TU[nc_TU$location == type]
    mtch <- findOverlaps(ChromHMM, tmp_TU)
    mtch.list <- split( queryHits(mtch), subjectHits(mtch) )
    
    for(i in as.numeric(names(mtch.list)))
    {
      idx <- sapply(mtch.list[[i]],
                    function(x) intersect(ranges(ChromHMM[x]), ranges(tmp_TU[i])) %>% width) %>% which.max
      tmp_state <- ChromHMM$name[mtch.list[[i]][idx]]
      nc_TU$ChromHMM[nc_TU$location == type][i] <- tmp_state
      chrm_TUs$Hits[ chrm_TUs$Type == type & chrm_TUs$State == tmp_state ] %+=% 1
    }
  }
  
  chrm_TUs$Type <- factor(chrm_TUs$Type, levels = c("intergenic","asRNA", "uaRNA", "conRNA", "daRNA", "usRNA"))
  chrm_TUs <- chrm_TUs[chrm_TUs$Type != "usRNA", ]
  chrm_TUs <- chrm_TUs[chrm_TUs$State %ni% c("12_LowSignal/RepetitiveElements","9_TranscriptionTransition", "1_Insulator"), ]
  chrm_TUs$State <- gsub(".*\\_", "\\1", chrm_TUs$State)
  chrm_TUs$State <- factor(chrm_TUs$State,
                           levels = c("WeakEnhancer", "Enhancer", "StrongEnhancer", "ActivePromoter", "TranscriptionElongation",
                                      "BivalentChromatin", "Intergenic", "RepressedChromatin", "Heterochromatin") )
  
  pdf("figs/ncTU ChromHMM hits.pdf", 6, 6)
  ggplot(chrm_TUs[chrm_TUs$Type %in% c("intergenic","asRNA"), ],
         aes(y = Hits, axis1 = Type, axis2 = State)) +
    geom_alluvium(aes(fill = Type), width = 1/14) +
    geom_stratum(width = 1/14, fill = colors_20[15], color = "black") +
    scale_x_discrete(limits = c("TU", "ChromHMM"), expand = c(.05, .1)) +
    geom_label(stat = "stratum", infer.label = TRUE, label.size = 0,
               nudge_x = c(0.02, 0.05, -0.10, -0.13, -0.04, -0.11, -0.15, -0.09, -0.10, -0.04, -0.09)) +
    scale_fill_manual(values = colors_20[c(2,5)]) +
    theme_setting +
    theme(panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")
  dev.off()
  
  # half of asRNAs associte with "promoter" or "enhancer"
  sum(chrm_TUs[chrm_TUs$Type == "asRNA", ][c(3,5,7,8), 3]) /
    sum(chrm_TUs[chrm_TUs$Type == "asRNA", 3])
}

