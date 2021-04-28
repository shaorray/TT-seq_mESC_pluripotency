setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# TTseq vs GROseq, PROseq
TU_GRO <- c(import.gff3("../data/TU_anno/other/TU_filter_2016_Flynn_GRO-seq_mES_SRR2036965.gtf"),
            import.gff3("../data/TU_anno/other/TU_filter_2017_Dorighi_GRO-seq_WT_SRR5466770.gtf"))
TU_PRO <- c(import.gff3("../data/TU_anno/other/TU_filter_2016_Jesse_PRO-Seq_SRR4041365.gtf"),
            import.gff3("../data/TU_anno/other/TU_filter_2018_Marta_ESC_PRO-Seq_SRR7300121.gtf"))
TU_TT <- readRDS("../fig1/data/TU.DE.mm10.RData")

TU_GRO_intergenic <- reduce(TU_GRO[TU_GRO$location == "intergenic"], min.gap = 1500)
TU_GRO_antisense <- reduce(TU_GRO[TU_GRO$location == "asRNA"], min.gap = 1500)

TU_PRO_intergenic <- reduce(TU_PRO[TU_PRO$location == "intergenic"], min.gap = 1500)
TU_PRO_antisense <- reduce(TU_PRO[TU_PRO$location == "asRNA"], min.gap = 1500)

TU_TT_intergenic <- reduce(TU_TT[TU_TT$location == "intergenic"], min.gap = 1500)
TU_TT_antisense <- reduce(TU_TT[TU_TT$location == "antisense"], min.gap = 1500)

TU_intergenic <- reduce(c(TU_GRO_intergenic, TU_PRO_intergenic, TU_TT_intergenic))
TU_antisense <- reduce(c(TU_GRO_antisense, TU_PRO_antisense, TU_TT_antisense))


TU_intergenic$Type[queryHits(findOverlaps(TU_intergenic, TU_GRO_intergenic))] <- "GRO-seq"
TU_intergenic$Type[queryHits(findOverlaps(TU_intergenic, TU_PRO_intergenic))] <- "PRO-seq"
TU_intergenic$Type[queryHits(findOverlaps(TU_intergenic, TU_TT_intergenic))] <- "TT-seq"
TU_intergenic$Type[queryHits(findOverlaps(TU_intergenic,
                                          intersect(TU_TT_intergenic, intersect(TU_GRO_intergenic, TU_PRO_intergenic))))] <- "Common"

TU_antisense$Type[queryHits(findOverlaps(TU_antisense, TU_GRO_antisense))] <- "GRO-seq"
TU_antisense$Type[queryHits(findOverlaps(TU_antisense, TU_PRO_antisense))] <- "PRO-seq"
TU_antisense$Type[queryHits(findOverlaps(TU_antisense, TU_TT_antisense))] <- "TT-seq"
TU_antisense$Type[queryHits(findOverlaps(TU_antisense,
                                          intersect(TU_TT_antisense, intersect(TU_GRO_antisense, TU_PRO_antisense))))] <- "Common"

TU_intergenic$GENCODE = TU_intergenic$FANTOM5 = TU_intergenic$NONCODE = TU_intergenic$ChromHMM = "None"
TU_antisense$GENCODE = TU_antisense$FANTOM5 = TU_antisense$NONCODE = TU_antisense$ChromHMM = "None"

# compare with known annotations -------------
# load reference annotations
gencode.gr = importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.vM25.annotation.gtf")
nc.gencode.gr = gencode.gr[gencode.gr$gene_type!='protein_coding' & gencode.gr$type=='gene']

noncode.gr = importRanges('../data/NONCODEv5_mm10.lncAndGene.bed')

enhancer.gr = importRanges('../data/F5.mm10.enhancers.bed')

ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name[ChromHMM$name %in% c("Enhancer", "StrongEnhancer", "WeakEnhancer")] <- "Enhancer"
ChromHMM$name[ChromHMM$name %in% c("TranscriptionElongation", "TranscriptionTransition")] <- "Transcription"


# TU overlapping with known references
# intergenic 
TU_intergenic$GENCODE[queryHits(findOverlaps(TU_intergenic, nc.gencode.gr))] <- "Known"
TU_intergenic$FANTOM5[queryHits(findOverlaps(TU_intergenic, enhancer.gr))] <- "Enhancer"
TU_intergenic$NONCODE[queryHits(findOverlaps(TU_intergenic, noncode.gr))] <- "Known"

mtch <- findOverlaps(TU_intergenic, ChromHMM)
TU_intergenic$ChromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]

# asRNA 
TU_antisense$GENCODE[queryHits(findOverlaps(TU_antisense, nc.gencode.gr))] <- "Known"
TU_antisense$FANTOM5[queryHits(findOverlaps(TU_antisense, enhancer.gr))] <- "Enhancer"
TU_antisense$NONCODE[queryHits(findOverlaps(TU_antisense, noncode.gr))] <- "Known"

mtch <- findOverlaps(TU_antisense, ChromHMM)
TU_antisense$ChromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]

# plot intergenic TUs overlap ----------------------------------------------------
library(Vennerable)
intergenic_hit <- cbind("TT-seq" = findOverlaps(TU_intergenic, TU_TT_intergenic) %>%
                          countQueryHits() > 0,
                        "PRO-seq" = findOverlaps(TU_intergenic, TU_PRO_intergenic) %>% 
                          countQueryHits() > 0,
                        "GRO-seq" = findOverlaps(TU_intergenic, TU_GRO_intergenic) %>%
                          countQueryHits() > 0)

sum(rowSums(intergenic_hit) == 3) # common 8001
sum(rowSums(intergenic_hit[, 1:2]) == 0) # GRO-seq 186191
sum(rowSums(intergenic_hit[, c(1, 3)]) == 0) # PRO-seq 9398
sum(rowSums(intergenic_hit[, 2:3]) == 0) # TT-seq 3804

sum(rowSums(intergenic_hit[, 1:2]) == 2 & !intergenic_hit[, 3]) # TT-seq & PRO-seq 280
sum(rowSums(intergenic_hit[, c(1, 3)]) == 2 & !intergenic_hit[, 2]) # TT-seq & GRO-seq 4572
sum(rowSums(intergenic_hit[, 2:3]) == 2 & !intergenic_hit[, 1]) # GRO-seq & PRO-seq 17960


# plot the fractions by alluvial ----------------------------------------------------
pacman::p_load(ggalluvial)

TU_intergenic_counts <- mcols(TU_intergenic) %>% as.data.frame() %>%
  dplyr::count(Type, GENCODE, FANTOM5, NONCODE, ChromHMM, sort = TRUE)

TU_intergenic_counts$Type <- factor(TU_intergenic_counts$Type, 
                                    levels = c("Common", "TT-seq", "GRO-seq", "PRO-seq"))
TU_intergenic_counts$ChromHMM <- factor(TU_intergenic_counts$ChromHMM, 
                                    levels = c("ActivePromoter", "Transcription", "BivalentChromatin",
                                               "Enhancer", "Insulator", "Heterochromatin", "Intergenic",
                                               "RepressedChromatin", "None"))

ggplot(TU_intergenic_counts,
       aes(y = n, axis1 = Type, axis2 = ChromHMM, axis3 = FANTOM5, axis4 = GENCODE, axis5 = NONCODE)) +
  geom_alluvium(aes(fill = Type), width = 1/14) +
  geom_stratum(width = 1/10, fill = colors_20[15], color = "black") +
  scale_x_discrete(limits = c("TU", "ChromHMM", "FANTOM5", "GENCODE", "NONCODE"), 
                   expand = c(.05, .1)) +
  geom_label(stat = "stratum", infer.label = TRUE, label.size = 0,
             size = 4, label.padding = unit(0.1, "lines"),) +
  scale_fill_manual(values = colors_9[c(8, 4:1)]) +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
ggsave(filename = "FigS1_intergenic_nascent_TUs_overlaps.png", 
       path = "../figS1/figs", width = 6, height = 6)

# antisense 
TU_antisense_counts <- mcols(TU_antisense) %>% as.data.frame() %>%
  dplyr::count(Type, GENCODE, FANTOM5, NONCODE, ChromHMM, sort = TRUE)

TU_antisense_counts$Type <- factor(TU_antisense_counts$Type, 
                                    levels = c("Common", "TT-seq", "GRO-seq", "PRO-seq"))
TU_antisense_counts$ChromHMM <- factor(TU_antisense_counts$ChromHMM, 
                                        levels = c("ActivePromoter", "Transcription", "BivalentChromatin",
                                                   "Enhancer", "Insulator", "Heterochromatin", "Intergenic",
                                                   "RepressedChromatin", "None"))

ggplot(TU_antisense_counts,
       aes(y = n, axis1 = Type, axis2 = ChromHMM, axis3 = FANTOM5, axis4 = GENCODE, axis5 = NONCODE)) +
  geom_alluvium(aes(fill = Type), width = 1/14) +
  geom_stratum(width = 1/10, fill = colors_20[16], color = "black") +
  scale_x_discrete(limits = c("TU", "ChromHMM", "FANTOM5", "GENCODE", "NONCODE"), 
                   expand = c(.05, .1)) +
  geom_label(stat = "stratum", infer.label = TRUE, label.size = 0,
             size = 4, label.padding = unit(0.1, "lines"),) +
  scale_fill_manual(values = colors_9[c(7, 4:1)]) +
  theme_setting +
  theme(panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
ggsave(filename = "FigS1_antisense_nascent_TUs_overlaps.png", 
       path = "../figS1/figs", width = 6, height = 6)

# load ATAC-seq peaks anno, Atlasi el al 2019 NCB, build on mm9 ----------------------------------------------------------------------------
ATAC.table = readxl::read_excel(path = "../data/41556_2019_310_MOESM6_ESM.xlsx", sheet = 1)
Dynamic_E1_E2_E3.table = readxl::read_excel(path = "../data/41556_2019_310_MOESM6_ESM.xlsx", sheet = 3)

ATAC.gr <- GRanges(seqnames = ATAC.table$chr, 
                   IRanges(start = ATAC.table$start, end = ATAC.table$end))
mcols(ATAC.gr) <- ATAC.table[, c("baseMeanA_FCS", "baseMeanB_2iDay1", "log2FoldChange_2iDay1", 
                                 "baseMeanB_2iDay7", "log2FoldChange_2iDay7")]

Dynamic_E1_E2_E3.gr <- GRanges(seqnames = Dynamic_E1_E2_E3.table$chr1, 
                               IRanges(start = Dynamic_E1_E2_E3.table$start, end = Dynamic_E1_E2_E3.table$end),
                               Type = Dynamic_E1_E2_E3.table$Type)

# load STARR-seq anno, Peng el al 2020 GB, build on mm9 ----------------------------------------------------------------------------
STARR.table = readxl::read_excel(path = "../data/13059_2020_2156_MOESM4_ESM.xlsx", sheet = 1)
STARR.gr <- GRanges(seqnames = STARR.table$chr, 
                    IRanges(start = STARR.table$start, end = STARR.table$end))
mcols(STARR.gr) <- STARR.table[, c("starr_peak_id", "enrichment_2iL", "enrichment_SL", "feature", "class", "GC_center_500bp")]

# lift to mm10
ATAC.mm10.gr <- liftOver(ATAC.gr, chain = import.chain("../data/liftOver_chains/mm9ToMm10.over.chain")) %>% unlist()
Dynamic_E1_E2_E3.mm10.gr <- liftOver(Dynamic_E1_E2_E3.gr, 
                                     chain = import.chain("../data/liftOver_chains/mm9ToMm10.over.chain")) %>% unlist()
STARR.mm10.gr <- liftOver(STARR.gr,
                          chain = import.chain("../data/liftOver_chains/mm9ToMm10.over.chain")) %>% unlist()

# save enhancer annotations
F5_enhancer <- importRanges("../data/F5.mm10.enhancers.bed")
F5_enhancer$name <- "FANTOM5"

super_enhancer <- importRanges("../data/super_enhancer_mm10.bed")
super_enhancer$name <- "Super_enhancer"

ChromHMM.mm10 <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM.mm10.enh <- ChromHMM.mm10[grep("Enhancer", ChromHMM.mm10$name)]
ChromHMM.mm10.enh$name <- paste0("ChromHMM_", gsub(".*_", "\\2", ChromHMM.mm10.enh$name))

STARR_enh.mm10 <- STARR.mm10.gr[STARR.mm10.gr$class %in% c("C1", "C2")] 
STARR_enh.mm10$name <- STARR_enh.mm10$class

# separate into 3 types: STARR with enhancer state / FANTOM5, STARR only, or enhancer state only
active_idx <- findOverlaps(STARR_enh.mm10, 
                           c(granges(F5_enhancer), granges(ChromHMM.mm10.enh))) %>%
              countQueryHits() > 0

enhancer_active_mm10 <- STARR_enh.mm10[active_idx]
enhancer_clean_mm10 <- STARR_enh.mm10[!active_idx]
enhancer_pseudo_mm10 <- ChromHMM.mm10.enh[findOverlaps(ChromHMM.mm10.enh, STARR_enh.mm10) %>%
                                            countQueryHits() == 0]

enhancer_mm10 <- c(granges(enhancer_active_mm10), granges(enhancer_clean_mm10), 
                   granges(enhancer_pseudo_mm10))

enhancer_mm10$type <- c(rep("active", length(enhancer_active_mm10)),
                        rep("clean", length(enhancer_clean_mm10)),
                        rep("pseudo", length(enhancer_pseudo_mm10)))

enhancer_mm10$name <- c(enhancer_active_mm10$name, enhancer_clean_mm10$name, enhancer_pseudo_mm10$name)
saveRDS(enhancer_mm10, "../data/enhancer_types.mm10.RData")

# compare enhancer overlap percentage ----------------------------------------------------------------------------
TU_intergenic$ATAC <- "None"
mtch <- findOverlaps(TU_intergenic, ATAC.mm10.gr)
TU_intergenic$ATAC[queryHits(mtch)] <- "ATAC_peak"
TU_intergenic$ATAC_baseMeanA_FCS <- NA
TU_intergenic$ATAC_baseMeanA_FCS[queryHits(mtch)] <- ATAC.mm10.gr$baseMeanA_FCS[subjectHits(mtch)]
TU_intergenic$ATAC_baseMeanA_2iDay1 <- NA
TU_intergenic$ATAC_baseMeanA_2iDay1[queryHits(mtch)] <- ATAC.mm10.gr$baseMeanB_2iDay1[subjectHits(mtch)]
TU_intergenic$ATAC_log2FoldChange_2iDay1 <- NA
TU_intergenic$ATAC_log2FoldChange_2iDay1[queryHits(mtch)] <- ATAC.mm10.gr$log2FoldChange_2iDay1[subjectHits(mtch)]
TU_intergenic$ATAC_log2FoldChange_2iDay7 <- NA
TU_intergenic$ATAC_log2FoldChange_2iDay7[queryHits(mtch)] <- ATAC.mm10.gr$log2FoldChange_2iDay7[subjectHits(mtch)]

TU_intergenic$Dynamic_E1_E2_E3 <- "None"
mtch <- findOverlaps(TU_intergenic, Dynamic_E1_E2_E3.mm10.gr)
TU_intergenic$Dynamic_E1_E2_E3[queryHits(mtch)] <- Dynamic_E1_E2_E3.mm10.gr$Type[subjectHits(mtch)]

TU_intergenic$STARR_class <- "None"
mtch <- findOverlaps(TU_intergenic, STARR.mm10.gr)
TU_intergenic$STARR_class[queryHits(mtch)] <- STARR.mm10.gr$class[subjectHits(mtch)]

# read Pol I/II/III density ---------------------------------------------------------------------------
TU_intergenic$Pol1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TX/mm10", 
                                          pattern = "Pol1", full.names = T), 
                               TU_intergenic + 50, fast = T) %>% rowMeans() %>% unlist() %>% unname()
TU_intergenic$Pol2 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TX/mm10", 
                                          pattern = "Pol2", full.names = T), 
                               TU_intergenic + 50, fast = T) %>% rowMeans() %>% unlist() %>% unname()
TU_intergenic$Pol3 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TX/mm10", 
                                          pattern = "Pol3", full.names = T), 
                               TU_intergenic + 50, fast = T) %>% rowMeans() %>% unlist() %>% unname()
TU_intergenic$Pol3 <- TU_intergenic$Pol3 / 5 # scale down since high background

Pol_mat <- mcols(TU_intergenic)[, c("Pol1", "Pol2", "Pol3")] %>% as.data.frame() 
Pol_mat[rowSums(Pol_mat < 0.05) == 3 | rowMax(Pol_mat) < 0.1, ] <- c(1, 1, 1)
Pol_mat <- Pol_mat / rowSums(Pol_mat) * 100
Pol_mat <- cbind(Pol_mat, 
                 ifelse(Pol_mat[, 1] > 90, 
                        "Pol1", ifelse(Pol_mat[, 2] > 90, 
                                       "Pol2", ifelse(Pol_mat[, 3] > 90,
                                                      "Pol3", ifelse(exp(rowSums(log(Pol_mat))) < 14000,
                                                                     "Mixed", "Low")))))
table(Pol_mat[, 4])
colnames(Pol_mat) <- c("Pol1", "Pol2", "Pol3", "Class")
Pol_mat$Class <- factor(Pol_mat$Class, c("Pol1", "Pol2", "Pol3", "Mixed", "Low"))

Pol_mat$Pol3_cofactors <- "False"
Pol3_idx <- sapply(list.files("../data/MACS2_peaks", "BRF1|TFIII", full.names = T), 
                   function(x) findOverlaps(TU_intergenic, importRanges(x)) %>% queryHits()) %>% unlist() %>% unique()
Pol_mat$Pol3_cofactors[Pol3_idx] <- "True"


# barplots ----------------------------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)

g1.1 <- data.frame(Overlap_rate = (table(TU_intergenic$Type, TU_intergenic$FANTOM5) / 
                                   c(table(TU_intergenic$Type)) * 100 )[, 1],
           Type = c("Common", "GRO-seq", "PRO-seq", "TT-seq")) %>%
  ggplot(aes(x = Type, y = Overlap_rate, fill = Type)) +
  geom_bar(stat="identity") +
  # scale_fill_manual(values = "grey75") +
  scale_y_continuous(breaks = c(0:5 * 10)) +
  coord_flip() +
  xlab("") + ylab("\nFANTOM5 enhancer ratio") +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.text = element_text(size = 11), 
        axis.title = element_text(size = 12))

g1.2 <- reshape2::melt((table(TU_intergenic$Type,
                      gsub("_.*", "\\2", TU_intergenic$Dynamic_E1_E2_E3)) /
                  c(table(TU_intergenic$Type)) * 100 )[, -4]) %>%
  ggplot(aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors_20[c(9,6,1)]) +
  coord_flip() +
  xlab("") + ylab("\nATAC-seq peaks ratio") + labs(fill = "Class") +
  theme_minimal() +
  theme(legend.position = "top", panel.grid = element_blank(),
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))

g1.3 <- reshape2::melt((table(TU_intergenic$Type, TU_intergenic$STARR_class) /
                  c(table(TU_intergenic$Type)) * 100 )[, -4]) %>%
  ggplot(aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_20[c(7,16,10)]) +
  coord_flip() +
  xlab("") + ylab("\nSTARR-seq peaks ratio") + labs(fill = "Class") +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))


g1.4 <- reshape2::melt((table(TU_intergenic$Type, Pol_mat$Class) /
                  c(table(TU_intergenic$Type)) * 100)) %>%
  ggplot(aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("brown1","dodgerblue1","limegreen", "grey50", "grey75")) +
  coord_flip() +
  xlab("") + ylab("\nPolymerase type enrichment") + labs(fill = "Class") +
  theme_minimal() +
  theme(legend.position = "top", panel.grid = element_blank(),
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12))

ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, g1.4, nrow=1, widths = c(6, 5, 5, 5)),
       filename = "FigS1_intergenic_nascent_RNA_classification.png",
       device = "png", path = "../figS1/figs", width = 17, height = 4)


Pol_den <- mcols(TU_intergenic)[, c("Pol1", "Pol2", "Pol3")] %>% as.data.frame()
Pol_den$Type <- TU_intergenic$Type
Pol_den <- reshape2::melt(Pol_den)
Pol_den$Type <- factor(Pol_den$Type, rev(c("Common", "GRO-seq", "PRO-seq", "TT-seq")))
Pol_den$value <- log(Pol_den$value)
Pol_den <- Pol_den[!is.infinite(Pol_den$value), ]
Pol_den$compare <- ifelse(Pol_den$Type == "Common", "Common", "Specific")
library(ggpubr)

compare_means(value ~ Type, group.by = "compare", data = Pol_den)

g2.1 <- ggplot(Pol_den, aes(x = Type, y = value, fill = Type)) +
  geom_violin(position = position_dodge(width = 0.8))+
  geom_boxplot(outlier.size = 0, width=.3, position = position_dodge(width = 0.8)) +
  facet_grid(.~variable) +
  stat_compare_means( aes(label = ..p.signif..), ref.group = "Common",
                      label.x = 2, label.y = 6) +
  xlab("\nTotal intergenic TUs") + ylab("Log density") +
  scale_fill_manual(values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D")) +
  ylim(c(-6, 6)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.title = element_text(size = 12), 
        strip.background = element_rect(fill="grey80")) 
ggsave(plot = g2.1,
       filename = "FigS1_intergenic_RNA_pol_density.png",
       device = "png", path = "../figS1/figs", width = 5, height = 4)

library(ggtern)
g2.2 <- ggtern(data = Pol_mat, 
               aes(x = Pol1, y = Pol2, z = Pol3, 
                   color = Class, 
                   shape = Pol3_cofactors),
       inherit.aes = FALSE) + 
  geom_mask() +
  geom_point(aes(size = c(0.1, 0.11)[as.numeric(as.factor(Pol3_cofactors))],
                 alpha = c(0.1, 0.5)[as.numeric(as.factor(Pol3_cofactors))])) +
  scale_color_manual(values = c("brown1","dodgerblue1","limegreen", "grey50", "grey75")) +
  scale_shape_manual(values = c(4, 1)) +
  theme_bw() +
  theme(legend.key.size = unit(0.6, "cm")) + 
  guides(color = guide_legend(override.aes = list(size = 2)),
         shape = guide_legend(override.aes = list(size = 2)),
         size = FALSE,
         alpha = FALSE)

ggsave(plot = g2.2,
       filename = "FigS1_intergenic_RNA_classification_triangle.png",
       device = "png", path = "../figS1/figs", width = 5, height = 4)

devtools::unload("ggtern")
R.methodsS3::setMethodS3("print", "ggplot", ggplot2:::print.ggplot)

# TTseq TU Pol2 type classification with GSE48895 ------------------------------------------------------------------------------------------
# filter 50min Trp sentitive TU -> Pol II TU 
# and 50min Flv sensitive TU -> Paused TU, mm9
TU.DE.mm9.gr <- readRDS("../fig1/data/TU.DE.mm9.RData")

# add features to annotated TUs
TU.DE.mm9.feature <- TU.DE.mm9.gr
mcols(TU.DE.mm9.feature) <- mcols(TU.DE.mm9.feature)[, 1:2]
TU.DE.mm9.feature$ATAC_peak <- findOverlaps(ATAC.gr, promoters(TU.DE.mm9.feature, upstream = 1000, downstream = 1000)) %>% 
                               countSubjectHits() > 0
TU.DE.mm9.feature$STARR_enh <- findOverlaps(Dynamic_E1_E2_E3.gr[!grepl("E3", Dynamic_E1_E2_E3.gr$Type)],
                                            promoters(TU.DE.mm9.feature, upstream = 1000, downstream = 1000)) %>% 
                               countSubjectHits() > 0
TU.DE.mm9.feature$ES_super_enh <- findOverlaps(importRanges("../data/ESC_SuperEnhancers.mm9.bed"),
                                               promoters(TU.DE.mm9.feature, upstream = 1000, downstream = 1000)) %>% 
                                  countSubjectHits() > 0

TU.DE.mm9.feature$Pol3 <- sapply(list.files("../data/MACS2_peaks/pol3_mm9", "RPC", full.names = T),
                                 function(x) import(x)) %>% Reduce(c, .) %>%
                          reduce() %>% findOverlaps(query = ., subject = TU.DE.mm9.feature) %>% 
                          countSubjectHits() > 0


# remove protein coding gene (0, 400) TSS region
TU.DE.mm9.gene_tss <- promoters(TU.DE.mm9.gr[TU.DE.mm9.gr$location == "protein_coding"],
                                upstream = 0, downstream = 400)
TU.DE.mm9.gene_body <- TU.DE.mm9.gr
TU.DE.mm9.gene_body[TU.DE.mm9.gr$location == "protein_coding"] <-
  foreach(i = seq_along(TU.DE.mm9.gene_tss), .combine = c) %dopar% {
  tmp = TU.DE.mm9.gene_tss[i]
  if (width(tmp) > 400) {
    if (as.character(strand(tmp)) == "+") {
      start(tmp) = end(TU.DE.mm9.gene_tss[i])
    } else {
      end(tmp) = start(TU.DE.mm9.gene_tss[i])
    }
  }
  tmp
}

dat <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2014_Jonker", ".bw$", full.names = T), 
                TU.DE.mm9.gene_body, fast = F) %>% 
        sweep(MARGIN =  2, c(-1, 1, -1, 1, -1, 1), "*")

dat_stranded <- dat[, c(1, 3, 5)]
dat_stranded[which(strand(TU.DE.mm9.gr) == "+"), ] <- dat[which(strand(TU.DE.mm9.gr) == "+"), c(2, 4, 6)]
dat_stranded[is.na(dat_stranded)] <- 0

# use the top 100 Pol III target genes as the references to normalise treatment conditions
idx_pol3 <- findOverlaps(TU.DE.mm9.gr, Pol3_peaks.mm9.gr) %>% queryHits()
idx_pol3 <- idx_pol3[as.integer(rank(rowSums(dat_stranded[idx_pol3, ]))) > (length(idx_pol3) - 100)] 

dat_stranded <- sweep(dat_stranded, MARGIN = 2, 
                      SizeFactorCal(dat_stranded[idx_pol3, ]),
                      FUN = "/")

# transform to log2 fold change
dat_cmp <- data.frame(Flv_change = log2(dat_stranded[, 1] / dat_stranded[, 3]),
                      Trp_change = log2(dat_stranded[, 2] / dat_stranded[, 3]))


plot_Trp_Flv_gro_seq <- function(dat_cmp, idx, show_y_text = F,
                                 Pol3_idx = NULL, cmp_idx = NULL, title = "NA") {
  
  dat_cmp_tmp <- dat_cmp[idx, ]
  idx_down_Flv <- dat_cmp_tmp[, 1] < (-2)
  idx_down_Trp <- dat_cmp_tmp[, 2] < (-2)
  
  rm_idx <- complete.cases(dat_cmp_tmp) & !is.infinite(rowSums(dat_cmp_tmp))
  dat_cmp_tmp <- dat_cmp_tmp[rm_idx, ]
  
  types <- rep(3, nrow(dat_cmp_tmp))
  types[which(idx_down_Trp[rm_idx])] <- 2
  types[which((idx_down_Flv & idx_down_Trp)[rm_idx])] <- 1
  
  down_Trp <- round(sum(idx_down_Trp & !idx_down_Flv, na.rm = T) / length(idx_down_Flv) * 100, 1)
  down_Trp_Flv <- round(sum(idx_down_Flv & idx_down_Trp, na.rm = T) / length(idx_down_Flv) * 100, 1)
  
  coord_rect <- data.frame(x1=c(0, 0, 0, 0),
                           x2=c(Inf, -Inf, Inf, -Inf),
                           y1=c(0, 0, 0, 0), 
                           y2=c(Inf, -Inf, -Inf, Inf), 
                           col=c('1','1','2','2'))
  
  g <- ggplot(dat_cmp_tmp, aes(x = Flv_change, y = Trp_change, 
                          alpha = with(dat_cmp_tmp, get_dens(Flv_change, Trp_change))^2,
                          color = types) ) +
    geom_rect(data = coord_rect, inherit.aes=FALSE,
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=col)) +
    geom_point(size = 0.5) +
    annotate("text", x = -4, y = -5.5, label = paste0(down_Trp_Flv, "%"), cex = 5) +
    annotate("text", x = 0, y = -5.5, label = paste0(down_Trp, "%"), cex = 5) +
    annotate("text", x = -6, y = 2, hjust = "left", vjust = "top",
             label = paste0("n = ", length(idx_down_Flv)), cex = 5) +
    scale_fill_manual(values=c('grey85','grey95')) +
    ggtitle(title) + xlab("") + ylab("") +
    xlim(c(-6, 3)) + ylim(c(-6, 2)) +
    theme_setting +
    theme(legend.position = "none",
          panel.border= element_blank(),
          axis.line = element_line(),
          plot.title = element_text(hjust = 0.5))
  
  if (show_y_text) {
    g <- g + ylab("GRO-seq log2FC Triptolide 50 min") 
  } else {
    g <- g + theme(axis.text.y = element_blank())
  }
  
  if (!is.null(Pol3_idx)) {
    g <- g + annotate("point", 
                      x = dat_cmp[intersect(which(idx), Pol3_idx), "Flv_change"],
                      y = dat_cmp[intersect(which(idx), Pol3_idx), "Trp_change"],
                      colour = "red", cex = 0.5)
  }
  
  if (!is.null(cmp_idx)) {
    g <- g + annotate("point", 
                      x = dat_cmp[intersect(which(idx), which(cmp_idx)), "Flv_change"],
                      y = dat_cmp[intersect(which(idx), which(cmp_idx)), "Trp_change"],
                      colour = "orange", cex = 0.3) + 
             annotate("text", x = 1, y = -4, hjust = "left",
                      label = paste("p <", signif(chisq.test(table(idx_down_Flv, cmp_idx[idx]))$p.value, 2)),
                      colour = "black", cex = 4)
  }
  g
}

g3.1 <- plot_Trp_Flv_gro_seq(dat_cmp = dat_cmp, idx = TU.DE.mm9.gr$location == "protein_coding", 
                     show_y_text = T, Pol3_idx = idx_pol3, cmp_idx = TU.DE.mm9.feature$ES_super_enh,
                     title = "mRNA")

g3.2 <- plot_Trp_Flv_gro_seq(dat_cmp = dat_cmp, idx = TU.DE.mm9.gr$location == "intergenic", 
                     show_y_text = F, Pol3_idx = idx_pol3, cmp_idx = TU.DE.mm9.feature$ES_super_enh,
                     title = "Intergenic")

g3.3 <- plot_Trp_Flv_gro_seq(dat_cmp = dat_cmp, idx = TU.DE.mm9.gr$location == "uaRNA", 
                     show_y_text = F, Pol3_idx = idx_pol3, cmp_idx = TU.DE.mm9.feature$ES_super_enh,
                     title = "uaRNA")

g3.4 <- plot_Trp_Flv_gro_seq(dat_cmp = dat_cmp, idx = TU.DE.mm9.gr$location == "antisense", 
                     show_y_text = F, Pol3_idx = idx_pol3, cmp_idx = TU.DE.mm9.feature$ES_super_enh,
                     title = "asRNA")

ggsave(plot = grid.arrange(g3.1, g3.2, g3.3, g3.4, nrow = 1, widths = c(5.1, 5, 5, 5)),
       filename = "FigS1_Trp_Flv_Gro-seq_change_classes.png",
       device = "png", path = "../figS1/figs", width = 15, height = 4)


# set thresholds at 2-fold changes
idx_Flv <- dat_stranded$GSE48895_V6.5_untreated / dat_stranded$GSE48895_V6.5_50minFP > 2
# table(idx)
table(idx_Flv, TU.DE.mm9.gr$location)

plot(trim_quantile(dat_stranded[TU.DE.mm9.gene_body$location == "protein_coding", c(1, 3)]),
     pch = 19, cex = 0.1)
abline(0, 1)
abline(0, 2)
points(dat_stranded[idx_Flv & TU.DE.mm9.gr$location == "protein_coding", c(1, 3)])


idx_Trp <- dat_stranded$GSE48895_V6.5_untreated / dat_stranded$GSE48895_V6.5_50minTrp > 2
table(idx_Trp, TU.DE.mm9.gr$location)

plot(trim_quantile(dat_stranded[TU.DE.mm9.gene_body$location == "protein_coding", c(2, 3)]), pch = 19, cex = 0.1)
abline(0, 1)
abline(0, 2)
points(dat_stranded[idx_Trp & TU.DE.mm9.gr$location == "protein_coding", c(2, 3)])





g2.3 <- ggplot(data.frame(value = table(TU_intergenic$STARR_class)),
               aes(x="", y=value.Freq, fill=value.Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  ylab("\nTotal STARR enhancers overlap") +
  scale_fill_manual(values = c(colors_20[c(7,16,10)], "grey75")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 0),
        axis.title = element_text(size = 0))