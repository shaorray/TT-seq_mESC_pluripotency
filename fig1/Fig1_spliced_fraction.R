# measure spliced reads contain in TTseq LRNA samples
# Rui Shao
# 2020 Feb

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
pacman::p_load(Rsubread, org.Mm.eg.db, TxDb.Mmusculus.UCSC.mm10.knownGene)

# system("mkdir figs")
# system("mkdir data")
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

# get genes intervals
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
gene.gr <- gene.gr[!duplicated(ranges(gene.gr))]

ks <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
res <- biomaRt::select(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL", 
                       columns = c("ENTREZID", "SYMBOL"))

gene.ann <- data.frame(GeneID = res$ENSEMBL[match(gene.gr$gene_id, res$ENTREZID)],
                       Chr = seqnames(gene.gr),
                       Start = start(gene.gr),
                       End = end(gene.gr),
                       Strand = strand(gene.gr))
gene.ann <- gene.ann[!is.na(gene.ann$GeneID), ]

# get sample sizes
spCounts <- readRDS("../data/spRPK_SL_2i.RData")
sp_L_mat <- spCounts[ ,grepl("LRNA",colnames(spCounts)) & ! grepl("PDS|SL2i",colnames(spCounts))]
colnames(sp_L_mat) <- gsub("(LRNA_)*","\\2",colnames(sp_L_mat))
LRNA.sizefactor <- SizeFactorCal(sp_L_mat)

# count total reads
bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10", pattern = ".bam$", full.names = T)
fc_PE <- Rsubread::featureCounts(bam_files, annot.ext=gene.ann, isPairedEnd=TRUE)
totalCounts <- setmm10LRNACounts(fc_PE$counts, LRNA.sizefactor)

# count spliced reads
spliced_bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10/spliced/", pattern = ".bam$", full.names = T)
fc_PE <- Rsubread::featureCounts(spliced_bam_files, annot.ext=gene.ann, isPairedEnd=TRUE)
splicedCounts <- setmm10LRNACounts(fc_PE$counts, LRNA.sizefactor)

# saveRDS(totalCounts, "data/totalCounts.RData")
# saveRDS(splicedCounts, "data/splicedCounts.RData")

# get spliced ratio
spliced_ratio <- splicedCounts / totalCounts
spliced_ratio <- spliced_ratio[rowSums(totalCounts) > 0, ] # remove inactive genes
spliced_ratio <- spliced_ratio[!apply(spliced_ratio, 1, function(x) any(is.na(x) | is.infinite(x))), ]
spliced_ratio <- spliced_ratio[apply(spliced_ratio, 1, function(x) sum(x > 0) > 4), ] # remove unspliced genes

colMedians(spliced_ratio)

cellCounts <- read.table("../data/cellCounts.txt", header = T)

# median spliced ratio with (0.25, 0.75) quantiles and SD
cellCounts$Spliced <- matrixStats::colMedians(spliced_ratio)[match(cellCounts$Samples, colnames(spliced_ratio))] 
cellCounts$Spliced_lower <- apply(spliced_ratio, 2, function(x) quantile(x, 0.25))[match(cellCounts$Samples, colnames(spliced_ratio))]
cellCounts$Spliced_upper <- apply(spliced_ratio, 2, function(x) quantile(x, 0.75))[match(cellCounts$Samples, colnames(spliced_ratio))]
cellCounts$Spliced_sd <- colSds(spliced_ratio)[match(cellCounts$Samples, colnames(spliced_ratio))]
cellCounts$Spliced_mean <-  colMeans(spliced_ratio)[match(cellCounts$Samples, colnames(spliced_ratio))] 

cellCounts$Color <- gsub("(.*)\\_rep.*", "\\1", cellCounts$Samples)

write.table(cellCounts, "data/Fig1_Splicing_cell_number.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

# statistical test
spliced_ratio2 <- spliced_ratio[, !grepl("SL_rep3", colnames(spliced_ratio))]
spliced_ratio_mean <- sapply(c("2i_2d", "mTORi_1d", "mTORi_2d", "SL"), 
                             function(x) rowMeans(spliced_ratio2[, grep(x, colnames(spliced_ratio))])) %>%
                      as.data.frame()
colMedians(as.matrix(spliced_ratio_mean))
colMeans(as.matrix(spliced_ratio_mean))

plot(density(log(spliced_ratio_mean$`2i_2d`/ spliced_ratio_mean$SL), na.rm = T))
plot(density(log(spliced_ratio_mean$mTORi_1d / spliced_ratio_mean$SL), na.rm = T))
plot(density(log(spliced_ratio_mean$mTORi_2d / spliced_ratio_mean$SL), na.rm = T))

median(inf.omit(na.omit((spliced_ratio_mean$`2i_2d` / spliced_ratio_mean$SL)))) # 0.9614648
median(inf.omit(na.omit((spliced_ratio_mean$mTORi_1d / spliced_ratio_mean$SL)))) # 0.7608696
median(inf.omit(na.omit((spliced_ratio_mean$mTORi_2d / spliced_ratio[, "SL_rep3"])))) # 1.100866

inf.omit = function(x) x[is.infinite(x)]
t.test(inf.omit(na.omit(log2(spliced_ratio_mean$SL / spliced_ratio_mean$`2i_2d`))))
t.test(inf.omit(na.omit(log2(spliced_ratio_mean$SL / spliced_ratio_mean$mTORi_1d))))
t.test(inf.omit(na.omit(log2(spliced_ratio[, "SL_rep3"] / spliced_ratio_mean$mTORi_2d))))

MASS::glm.nb(SL ~ `2i_2d`, data = spliced_ratio_mean) %>% summary()
MASS::glm.nb(SL ~ mTORi_1d, data = spliced_ratio_mean) %>% summary()
MASS::glm.nb(SL ~ mTORi_2d, data = spliced_ratio_mean) %>% summary()



# --------------------------------------------------------------------------------------------
cellCounts <- read.table("data/Fig1_Splicing_cell_number.txt", header = T)
cellCounts$Color <- factor(cellCounts$Color, c("SL", "2i_2d", "mTORi_1d", "mTORi_2d"))
cellCounts$Samples <- factor(cellCounts$Samples, unique(cellCounts$Samples))

# plot
ggplot(cellCounts, aes(x = Samples, y = Spliced * 100, fill = Color)) + 
  geom_bar(stat = "identity") +
  # geom_point(cex = 3, pch = 15) + 
  geom_linerange(aes(ymin=Spliced_lower * 100, ymax=Spliced_upper * 100)) +
  # xlab("\nCell density (million / cm2)") +
  xlab("") +
  ylab("Spliced %\n") +
  ylim(0, 4) +
  scale_fill_manual(values = colors_20[c(13, 2, 20, 7)]) +
  labs(fill='Samples') +
  theme_setting +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust=1))

ggsave(filename = "Fig1.Spliced_ratio_sample.pdf", path = "figs", device = "pdf",
       width = 5, height = 4)
# ggsave(filename = "Fig1.Splicing_cell_number.pdf", path = "figs", device = "pdf",
#        width = 5, height = 4)
