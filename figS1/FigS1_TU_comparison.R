setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# XCI dataset
TU_TT_1_0 <- importRanges("../data/TU_anno/other/TU_filter+2018_Jan_TTseq_TX1072_WT_0h.gtf")
TU_TT_1_4 <- importRanges("../data/TU_anno/other/TU_filter+2018_Jan_TTseq_TX1072_WT_4h.gtf")
TU_TT_1_8 <- importRanges("../data/TU_anno/other/TU_filter+2018_Jan_TTseq_TX1072_WT_8h.gtf")
TU_TT_1_12 <- importRanges("../data/TU_anno/other/TU_filter+2018_Jan_TTseq_TX1072_WT_12h.gtf")
TU_TT_1_24 <- importRanges("../data/TU_anno/other/TU_filter+2018_Jan_TTseq_TX1072_WT_24h.gtf")

# other nascent RNA datasets
TU_GRO_1 <- import.gff3("../data/TU_anno/TU_filter_2016_Flynn_GRO-seq_mES_SRR2036965.gtf")
TU_GRO_2 <- import.gff3("../data/TU_anno/TU_filter_2017_Dorighi_GRO-seq_WT_SRR5466770.gtf")
TU_PRO_1 <- import.gff3("../data/TU_anno/TU_filter_2016_Jesse_PRO-Seq_SRR4041365.gtf")
TU_PRO_2 <- import.gff3("../data/TU_anno/TU_filter_2018_Marta_ESC_PRO-Seq_SRR7300121.gtf")

# functions
jacIntersect <- function(tu1, tu2)
{ # compute jaccard correlation between two input TUs
  sum(width(intersect(tu1, tu2))) / sum(width(reduce(c(tu1, tu2))))
}


# novel transcripts
num_TUs <- rbind( cbind(TU = table(TU_TT_1_0$location), 
                        Len = aggregate( width(TU_TT_1_0) / 1000, list(TU_TT_1_0$location), sum)[,2],
                        Type = names(table(TU_TT_1_0$location)),
                        Sample = rep("XCI_0h", length(table(TU_TT_1_0$location)))),
                  
                  cbind(TU = table(TU_TT_1_4$location), 
                        Len = aggregate( width(TU_TT_1_4) / 1000, list(TU_TT_1_4$location), sum)[,2],
                        Type = names(table(TU_TT_1_4$location)),
                        Sample = rep("XCI_4h", length(table(TU_TT_1_4$location)))),
                  
                  cbind(TU = table(TU_TT_1_8$location), 
                        Len = aggregate( width(TU_TT_1_8) / 1000, list(TU_TT_1_8$location), sum)[,2],
                        Type = names(table(TU_TT_1_8$location)),
                        Sample = rep("XCI_8h", length(table(TU_TT_1_8$location)))),
                  
                  cbind(TU = table(TU_TT_1_12$location), 
                        Len = aggregate( width(TU_TT_1_12) / 1000, list(TU_TT_1_12$location), sum)[,2],
                        Type = names(table(TU_TT_1_12$location)),
                        Sample = rep("XCI_12h", length(table(TU_TT_1_12$location)))),
                  
                  cbind(TU = table(TU_TT_1_24$location), 
                        Len = aggregate( width(TU_TT_1_24) / 1000, list(TU_TT_1_24$location), sum)[,2],
                        Type = names(table(TU_TT_1_24$location)),
                        Sample = rep("XCI_24h", length(table(TU_TT_1_24$location)))),
                  
                  cbind(TU = table(TU_SL_LRNA$location), 
                        Len = aggregate( width(TU_SL_LRNA) / 1000, list(TU_SL_LRNA$location), sum)[,2],
                        Type = names(table(TU_SL_LRNA$location)),
                        Sample = rep("SL", length(table(TU_SL_LRNA$location)))),
                  
                  cbind(TU = table(TU_2i_LRNA$location), 
                        Len = aggregate( width(TU_2i_LRNA) / 1000, list(TU_2i_LRNA$location), sum)[,2],
                        Type = names(table(TU_2i_LRNA$location)),
                        Sample = rep("2i_2d", length(table(TU_2i_LRNA$location)))),
                  
                  cbind(TU = table(TU_mTORi_1d_LRNA$location), 
                        Len = aggregate( width(TU_mTORi_1d_LRNA) / 1000, list(TU_mTORi_1d_LRNA$location), sum)[,2],
                        Type = names(table(TU_mTORi_1d_LRNA$location)),
                        Sample = rep("mTORi_1d", length(table(TU_mTORi_1d_LRNA$location)))),
                  
                  cbind(TU = table(TU_mTORi_2d_LRNA$location), 
                        Len = aggregate( width(TU_mTORi_2d_LRNA) / 1000, list(TU_mTORi_2d_LRNA$location), sum)[,2],
                        Type = names(table(TU_mTORi_2d_LRNA$location)),
                        Sample = rep("mTORi_2d", length(table(TU_mTORi_2d_LRNA$location))))
) 
num_TUs[num_TUs[, 3] == "protein_coding", 3] <- "mRNA"
num_TUs[num_TUs[, 3] == "antisense", 3] <- "asRNA"
num_TUs <- num_TUs[!num_TUs[, 3] %in% c("sense_intragenic", "dsRNA"), ] # keep high confidence TU
num_TUs <- num_TUs %>% as.data.frame
num_TUs$TU <- num_TUs$TU %>% as.character %>% as.numeric
num_TUs$Len <- num_TUs$Len %>% as.character %>% as.numeric
num_TUs$Type <- factor(num_TUs$Type, c("intergenic", "mRNA", "asRNA", "uaRNA", "conRNA", "daRNA", "usRNA"))
num_TUs$Sample <- factor(num_TUs$Sample, c("XCI_0h", "XCI_4h", "XCI_8h", "XCI_12h", "XCI_24h", "SL", "2i_2d", "mTORi_1d", "mTORi_2d"))


pdf("figs/XCI_TU number by sample.pdf", 6, 5)
ggplot(num_TUs, aes(x = Type, y = TU, fill = Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("TU number") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(5, "Paired"), colors_20[c(13, 2, 20, 7)])) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

pdf("figs/XCI_TU total length.pdf", 6, 5)
ggplot(num_TUs[num_TUs$Type != "mRNA", ], aes(x = Type, y = Len, fill = Sample)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("TU total length (Kb)") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(5, "Paired"), colors_20[c(13, 2, 20, 7)])) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

# TPM
# -------
tpm_TUs <- rbind( cbind( TPM = as.numeric(TU_TT_1_0$expr) / sum(as.numeric(TU_TT_1_0$expr)) * 1e6,
                         Type = TU_TT_1_0$location,
                         Sample = rep("XCI_0h", length(TU_TT_1_0)) ),
                  
                  cbind( TPM = as.numeric(TU_TT_1_4$expr) / sum(as.numeric(TU_TT_1_4$expr)) * 1e6,
                         Type = TU_TT_1_4$location,
                         Sample = rep("XCI_4h", length(TU_TT_1_4)) ),
                  
                  cbind( TPM = as.numeric(TU_TT_1_8$expr) / sum(as.numeric(TU_TT_1_8$expr)) * 1e6,
                         Type = TU_TT_1_8$location,
                         Sample = rep("XCI_8h", length(TU_TT_1_8)) ),
                  
                  cbind( TPM = as.numeric(TU_TT_1_12$expr) / sum(as.numeric(TU_TT_1_12$expr)) * 1e6,
                         Type = TU_TT_1_12$location,
                         Sample = rep("XCI_12h", length(TU_TT_1_12)) ),
                  
                  cbind( TPM = as.numeric(TU_TT_1_24$expr) / sum(as.numeric(TU_TT_1_24$expr)) * 1e6,
                         Type = TU_TT_1_24$location,
                         Sample = rep("XCI_24h", length(TU_TT_1_24)) )
) 

tpm_TUs[tpm_TUs[, 2] == "protein_coding", 2] <- "mRNA"
tpm_TUs[tpm_TUs[, 2] == "antisense", 2] <- "asRNA"
tpm_TUs <- tpm_TUs[!tpm_TUs[, 2] %in% c("sense_intragenic", "dsRNA"), ] # keep high confidence TUs

tpm_TUs <- tpm_TUs %>% as.data.frame
tpm_TUs$Type <- factor(tpm_TUs$Type, c("mRNA", "intergenic","asRNA", "uaRNA", "conRNA", "daRNA", "usRNA"))
tpm_TUs$TPM <- tpm_TUs$TPM %>% as.character %>% as.numeric
tpm_TUs$Sample <- factor(tpm_TUs$Sample, c("XCI_0h", "XCI_4h", "XCI_8h", "XCI_12h", "XCI_24h"))

pdf("figs/XCI_TU_TPM.pdf", 9, 5)
ggplot(tpm_TUs, aes(x = Type, y = TPM , fill = Sample)) +
  geom_boxplot(outlier.size = 0) +
  xlab("") +
  ylab("TPM") +
  scale_y_log10() +
  ylim(c(0, 30)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Paired")) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

# TTseq vs GROseq, PROseq

# -------
tpm_TUs <- rbind( cbind( TPM = as.numeric(TU_TT_1_0$expr) / sum(as.numeric(TU_TT_1_0$expr)) * 1e6,
                         Type = TU_TT_1_0$location,
                         Sample = rep("TT-seq1", length(TU_TT_1_0)) ),
                  
                  cbind( TPM = as.numeric(TU_SL_LRNA$expr) / sum(as.numeric(TU_SL_LRNA$expr)) * 1e6,
                         Type = TU_SL_LRNA$location,
                         Sample = rep("TT-seq2", length(TU_SL_LRNA)) ),
                  
                  cbind( TPM = as.numeric(TU_PRO_1$expr) / sum(as.numeric(TU_PRO_1$expr)) * 1e6,
                         Type = TU_PRO_1$location,
                         Sample = rep("PRO-seq1", length(TU_PRO_1)) ),
                  
                  cbind( TPM = as.numeric(TU_PRO_2$expr) / sum(as.numeric(TU_PRO_2$expr)) * 1e6,
                         Type = TU_PRO_2$location,
                         Sample = rep("PRO-seq2", length(TU_PRO_2)) ),
                  
                  cbind( TPM = as.numeric(TU_GRO_1$expr) / sum(as.numeric(TU_GRO_1$expr)) * 1e6,
                         Type = TU_GRO_1$location,
                         Sample = rep("GRO-seq1", length(TU_GRO_1)) ),
                  
                  cbind( TPM = as.numeric(TU_GRO_2$expr) / sum(as.numeric(TU_GRO_2$expr)) * 1e6,
                         Type = TU_GRO_2$location,
                         Sample = rep("GRO-seq2", length(TU_GRO_2)) )
) 

tpm_TUs[tpm_TUs[, 2] == "protein_coding", 2] <- "mRNA"
tpm_TUs[tpm_TUs[, 2] == "antisense", 2] <- "asRNA"
tpm_TUs <- tpm_TUs[!tpm_TUs[, 2] %in% c("sense_intragenic", "dsRNA"), ] # keep high confidence TUs

tpm_TUs <- tpm_TUs %>% as.data.frame
tpm_TUs$Type <- factor(tpm_TUs$Type, c("mRNA", "intergenic","asRNA", "uaRNA", "conRNA", "daRNA", "usRNA"))
tpm_TUs$TPM <- tpm_TUs$TPM %>% as.character %>% as.numeric
tpm_TUs$Sample <- factor(tpm_TUs$Sample, c("TT-seq1", "TT-seq2", "PRO-seq1", "PRO-seq2", "GRO-seq1", "GRO-seq2"))

pdf("figs/Comparison_TU_TPM.pdf", 9, 5)
ggplot(tpm_TUs, aes(x = Type, y = TPM , fill = Sample)) +
  geom_boxplot(outlier.size = 0) +
  xlab("") +
  ylab("TPM") +
  scale_y_log10() +
  ylim(c(0, 30)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Paired")) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

# ---------

# reads distribution
rc_TUs <- rbind( cbind(aggregate(as.numeric(TU_TT_1_0$expr) * width(TU_TT_1_0), list(TU_TT_1_0$location), sum),
                       Sample = rep("TT-seq1", length(table(TU_TT_1_0$location)))),
                 cbind(aggregate(as.numeric(TU_SL_LRNA$expr) * width(TU_SL_LRNA), list(TU_SL_LRNA$location), sum),
                       Sample = rep("TT-seq2", length(table(TU_SL_LRNA$location)))),
                 cbind(aggregate(as.numeric(TU_PRO_1$expr) * width(TU_PRO_1), list(TU_PRO_1$location), sum),
                       Sample = rep("PRO-seq1", length(table(TU_PRO_1$location)))),
                 cbind(aggregate(as.numeric(TU_PRO_2$expr) * width(TU_PRO_2), list(TU_PRO_2$location), sum),
                       Sample = rep("PRO-seq2", length(table(TU_PRO_2$location)))),
                 cbind(aggregate(as.numeric(TU_GRO_1$expr) * width(TU_GRO_1), list(TU_GRO_1$location), sum),
                       Sample = rep("GRO-seq1", length(table(TU_GRO_1$location)))),
                 cbind(aggregate(as.numeric(TU_GRO_2$expr) * width(TU_GRO_2), list(TU_GRO_2$location), sum),
                       Sample = rep("GRO-seq2", length(table(TU_GRO_2$location))))
) 

rc_TUs[rc_TUs[, 1] == "protein_coding", 1] <- "mRNA"
rc_TUs[rc_TUs[, 1] == "antisense", 1] <- "asRNA"
rc_TUs <- rc_TUs[!rc_TUs[, 1] %in% c("sense_intragenic", "dsRNA"), ] # keep high confidence TUs
rc_TUs[rc_TUs[, 3] == "TT-seq1", 2] <- rc_TUs[rc_TUs[, 3] == "TT-seq1", 2] / sum(rc_TUs[rc_TUs[, 3] == "TT-seq1", 2]) * 100
rc_TUs[rc_TUs[, 3] == "TT-seq2", 2] <- rc_TUs[rc_TUs[, 3] == "TT-seq2", 2] / sum(rc_TUs[rc_TUs[, 3] == "TT-seq2", 2]) * 100
rc_TUs[rc_TUs[, 3] == "PRO-seq1", 2] <- rc_TUs[rc_TUs[, 3] == "PRO-seq1", 2] / sum(rc_TUs[rc_TUs[, 3] == "PRO-seq1", 2]) * 100
rc_TUs[rc_TUs[, 3] == "PRO-seq2", 2] <- rc_TUs[rc_TUs[, 3] == "PRO-seq2", 2] / sum(rc_TUs[rc_TUs[, 3] == "PRO-seq2", 2]) * 100
rc_TUs[rc_TUs[, 3] == "GRO-seq1", 2] <- rc_TUs[rc_TUs[, 3] == "GRO-seq1", 2] / sum(rc_TUs[rc_TUs[, 3] == "GRO-seq1", 2]) * 100
rc_TUs[rc_TUs[, 3] == "GRO-seq2", 2] <- rc_TUs[rc_TUs[, 3] == "GRO-seq2", 2] / sum(rc_TUs[rc_TUs[, 3] == "GRO-seq2", 2]) * 100
rc_TUs <- rc_TUs[!rc_TUs[, 1] %in% c("mRNA"), ] # keep non-coding TUs

colnames(rc_TUs) <- c("Type", "Abundance", "Sample")
rc_TUs <- rc_TUs %>% as.data.frame
# rc_TUs$TU <- rc_TUs$TU %>% as.character %>% as.numeric
rc_TUs$Type <- factor(rc_TUs$Type, c("intergenic", "asRNA", "uaRNA", "conRNA", "daRNA", "usRNA"))


pdf("figs/Comparison_TU_reads_abundance.pdf", 6, 5)
ggplot(rc_TUs, aes(x = Type, y = Abundance, fill = Sample)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.6) +
  xlab("") +
  ylab("Reads Abundance %") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Paired")) +
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

# compare with known intergenic TUs -------------
# load reference annotations
gencode.gr = importRanges('~/Documents/genomeDir/Mus_musculus.GRCm38.79.gtf')
nc.gencode.gr = gencode.gr[gencode.gr$gene_biotype!='protein_coding' & gencode.gr$type=='gene']

noncode.gr = importRanges('../data/NONCODEv5_mm10.lncAndGene.bed')

enhancer.gr = importRanges('../data/F5.mm10.enhancers.bed')

annoTU <- function(TU, ref.list, .names)
{
  match_vec <- rep("Unkown", length(TU))
  for (i in seq_along(ref.list)) 
  {
    idx <- findOverlaps(TU, ref.list[[i]]) %>% queryHits 
    match_vec[idx] <- .names[i]
  }
  return( table(match_vec) / length(TU))
}

Target_type = "intergenic"
Target_type = "asRNA"

ref_mtch <- rbind(cbind(Freq = annoTU(TU_TT_1_0[TU_TT_1_0$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("TT-seq1", 4)),
                  
                  cbind(Freq = annoTU(TU_SL_LRNA[TU_SL_LRNA$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("TT-seq2", 4)),
                  
                  cbind(Freq = annoTU(TU_PRO_1[TU_PRO_1$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("PRO-seq1", 4)),
                  
                  cbind(Freq = annoTU(TU_PRO_2[TU_PRO_2$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("PRO-seq2", 4)),
                  
                  cbind(Freq = annoTU(TU_GRO_1[TU_GRO_1$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("GRO-seq1", 4)),
                  
                  cbind(Freq = annoTU(TU_GRO_2[TU_GRO_2$location == Target_type], 
                                      list(nc.gencode.gr, noncode.gr, enhancer.gr), 
                                      c("GENCODE", "NONCODE", "Enhancer")),
                        Type = c("Enhancer", "GENCODE", "NONCODE", "Unkown"),
                        Sample = rep("GRO-seq2", 4))
)

ref_mtch <- ref_mtch %>% as.data.frame
ref_mtch$Freq <- ref_mtch$Freq %>% as.character %>% as.numeric %>% "*"(100)
ref_mtch$Sample <- factor(ref_mtch$Sample, c("TT-seq1", "TT-seq2", "PRO-seq1", "PRO-seq2", "GRO-seq1", "GRO-seq2"))


pdf("figs/Comparison_TU_intergenic_type.pdf", 6, 5)
ggplot(ref_mtch, aes(x = Sample, y = Freq, fill = Type)) +
  geom_bar(stat="identity", position=position_stack()) +
  xlab("") +
  ylab("Overlapping % \n") +
  scale_fill_manual(values = colors_9) +
  labs(title = "Intergenic TUs") + 
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()

pdf("figs/Comparison_TU_asRNA_type.pdf", 6, 5)
ggplot(ref_mtch, aes(x = Sample, y = Freq, fill = Type)) +
  geom_bar(stat="identity", position=position_stack()) +
  xlab("") +
  ylab("Overlapping % \n") +
  scale_fill_manual(values = colors_9) +
  labs(title = "asRNA TUs") + 
  theme_setting +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line())
dev.off()