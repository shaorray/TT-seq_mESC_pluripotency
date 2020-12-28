#-----------------------------------------------------------------------------------------------------
# Multiple genomic and epigenomic features analysis with tx, half-life and RNA copy number
#-----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# prepare TU intervals --------------------------------------------------------------------------------------------
if (T) {
  TU.DE.mm9.gr <- readRDS("../fig1/data/TU.DE.mm9.RData")
  TU.DE.mm9.gr <- TU.DE.mm9.gr[TU.DE.mm9.gr$baseMean_LRNA > 0]
  
  # clean up mRNA annotation (84 TUs removed)
  non_mRNA.idx <- (!is.na(TU.DE.mm9.gr$gene_id)) %>% which() %>%
    "["(biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                        keys = gsub("\\..*", "", TU.DE.mm9.gr$gene_id[!is.na(TU.DE.mm9.gr$gene_id)]),
                        keytype = "GENEID",
                        columns = "GENEBIOTYPE")$"GENEBIOTYPE" != "protein_coding") 
  
  # TSS intervals
  TU.TSS.mm9.gr <- promoters(TU.DE.mm9.gr, upstream = 500, downstream = 500)
  mcols(TU.TSS.mm9.gr) <- NULL
  H3K4me3.peaks <- import("../data/H3K4m3_SL_CTR_peaks.broadPeak")
  mcols(H3K4me3.peaks) <- NULL
  mtch <- findOverlaps(TU.TSS.mm9.gr, H3K4me3.peaks)
  TU.TSS.mm9.gr[queryHits(mtch)] <- H3K4me3.peaks[subjectHits(mtch)]
  
  # get ChromHMM states
  ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
  ChromHMM$name <- gsub("(.*?)_(.*)", "\\2", ChromHMM$name)
  ChromHMM <- liftOver(ChromHMM, import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>% unlist()
  mtch <- findOverlaps(promoters(TU.DE.mm9.gr, upstream = 500, downstream = 500),
                       ChromHMM)
  TU.TSS.mm9.gr$chromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]
  
  # get directionality
  TU.as <- promoters(TU.DE.mm9.gr, upstream = 2000, downstream = 500)
  levels(strand(TU.as)) <- c("-", "+", "*")
  mtch_dir <- findOverlaps(TU.DE.mm9.gr, TU.as)
  TU.TSS.mm9.gr$Direction <- ifelse(mtch_dir %>% countQueryHits() > 0,
                                    "Bidirection", "Unidirection" )
}

# match with tx parameter estimation
tu_names <- paste(TU.DE.mm9.gr$location, seqnames(TU.DE.mm9.gr), start(TU.DE.mm9.gr),
                  end(TU.DE.mm9.gr), strand(TU.DE.mm9.gr), sep = "_")
TU.DE.mm9.gr$gene_id <- gsub("\\..*", "", TU.DE.mm9.gr$gene_id)
tu_names[TU.DE.mm9.gr$location == "protein_coding"] <- TU.DE.mm9.gr$gene_id[TU.DE.mm9.gr$location == "protein_coding"]

sl_tu_counts_Rates <-
  readRDS("../figS2/data/sample_tx_tu_rates_combined.mm9.RData") %>% 
  dplyr::filter(Sample == "SL")

# process epigenome features of TUs --------------------------------------------------------------------------------------------
if (T) {
  TU_feature_table <- sl_tu_counts_Rates[match(tu_names, sl_tu_counts_Rates$gene_id), ]
  TU_feature_table$TX_RPK <- TU.DE.mm9.gr$baseMean_LRNA / width(TU.DE.mm9.gr) * 1e3
  
  # DNA elements on TSS --------------------------------------------------------------------------------------------
  TU_feature_table$CpG_num <- 
    interval_pattern_hits(intervals = promoters(TU.DE.mm9.gr, 
                                                upstream = 1000,
                                                downstream = 50),
                          pattern = "CG", 
                          which_genome = "mm9")
  
  TU_feature_table$TATA_num <-
    interval_pattern_hits(intervals = promoters(TU.DE.mm9.gr, 
                                                upstream = 1000,
                                                downstream = 50),
                          pattern = "TATA", # using c("TATAAAA", "TATATAA", "TATAAAT", "TATATAT") gives the results
                          which_genome = "mm9")
  
  # DNA methylation -----------------------------------------------------------------------------------------------------
  TU_feature_table$TSS_DNAme. <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/CpG/GSM1127953_DNA_CpG_E14_serum_rep.bw", 
                                   TU.TSS.mm9.gr, fast = T) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  TU_feature_table$GB_DNAme. <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/CpG/GSM1127953_DNA_CpG_E14_serum_rep.bw", 
                                           TU.DE.mm9.gr, fast = T) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  # Chromatin accessibility --------------------------------------------------------------------------------------------
  TU_feature_table$DHS <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/DHS", full.names = T), 
                                   TU.TSS.mm9.gr, fast = T) %>% rowMeans() %>% 
    unlist() %>% unname()
  TU_feature_table$ATAC <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/ATAC", "Dorighi_ATAC", full.names = T), 
                                    TU.TSS.mm9.gr, fast = T) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  # General transcription factors / chromatin remodelers on TSS -----------------------------------------------------------------------
  TU_feature_table$TBP <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2015_Langer_TBP_ESwt.mm9.bw", 
                                   TU.TSS.mm9.gr, fast = T) %>% 
    unlist() %>% unname()
  
  TU_feature_table$TFIIB <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2015_Langer_TFIIB_ESwt.mm9.bw", 
                                     TU.TSS.mm9.gr, fast = T) %>% 
    unlist() %>% unname()
  
  TU_feature_table$Cebpb <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2015_Goode_MES_Cebpb.mm9.bw", 
                                     TU.TSS.mm9.gr, fast = T) %>%
    unlist() %>% unname()
  
  TU_feature_table$Med1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2018_Sabari_MED1_Control.mm9.bw", 
                                    TU.TSS.mm9.gr, fast = T) %>% 
    unlist() %>% unname()
  
  TU_feature_table$Brd4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2018_Sabari_BRD4_Control.mm9.bw", 
                                    TU.TSS.mm9.gr, fast = T) %>% 
    unlist() %>% unname()
  
  TU_feature_table$CTCF <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/GSE92407", 
                                               "CTCF-ChIP-serum", full.names = T), 
                                    TU.TSS.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Med1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                               "Xue_IP_Med1", full.names = T), 
                                    TU.TSS.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$TFIID <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2012_Ku_TFIID.mm9.bw",
                                     TU.TSS.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$Sp1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                               "Sp1_rep", full.names = T), 
                                    TU.TSS.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Ezh2 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                              "EZH2", full.names = T), 
                                    TU.TSS.mm9.gr, fast = T) %>% rowMeans()
  
  # General transcription factor / chromatin reader / regulator on gene body -----------------------------------------------------------------------
  TU_feature_table$LaminB <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/domain/2017_Poleshko_ESC_LaminB_ChIP-seq.mm9.bw",
                                      TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$HP1a <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/domain", "HP1a", full.names = T),
                                    TU.DE.mm9.gr, fast = T) %>% rowMeans() %>% unname()
  
  TU_feature_table$Suz12 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                                "Kaneko_E14.SUZ12", full.names = T),
                                     TU.DE.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Aff4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2011_Lin_AFF4_0hr_ChIPSeq.mm9.bw",
                                    TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$Yy1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/GSE92407", 
                                              "Yy1-ChIP-serum", full.names = T), 
                                   TU.DE.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Chd1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                                "Chd1", full.names = T),
                                     TU.DE.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Chd2 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                               "Chd2", full.names = T),
                                    TU.DE.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Chd9 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                               "Chd9", full.names = T),
                                    TU.DE.mm9.gr, fast = T) %>% rowMeans()
  
  TU_feature_table$Ring1b <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/domain/GSM2393579_ES_Ring1b.density.bigWig",
                                    TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$Phc1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/domain/GSM2393580_ES_Phc1.density.bigWig",
                                      TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()

    
  # Pluripotent transcription factors on TSS --------------------------------------------------------------------------------------------
  tmp <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", pattern = "2016_Chronis_ESC", full.names = T), 
                  TU.TSS.mm9.gr, fast = T) #%>% '*'(width(TU.TSS.mm9.gr) / 1e3)
  TU_feature_table <- cbind(TU_feature_table,
                            `colnames<-`(tmp, gsub(".*ESC_(.*)_ChIP.*", "\\1", colnames(tmp))) )
  
  # Histone marks on gene bodies --------------------------------------------------------------------------------------------
  tmp <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/Histone/", ".bw$", full.names = T), 
                  TU.DE.mm9.gr, fast = T) 
  tmp <- tmp[, -19]
  colnames(tmp) <- gsub("_ChIP.*|_rep.*|_wt.*|_WT", "\\1", gsub(".*20.*_(H.*).mm9.bw", "\\1", colnames(tmp)))
  TU_feature_table <- cbind(TU_feature_table, tmp)
  
  TU_feature_table$H33_YFP <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/GSM2582412_ESC_H3.3WT_YFP.bigwig",
                                       TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  TU_feature_table$H2AZ_WT <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/GSM984545_H2AZ_WT.BigWig",
                                       TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$H2A <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/new/2020_Wen_H2a_mnase32_x_chip_rep1.mm9.bw",
                                   TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  TU_feature_table$H33 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/new/2020_Wen_H33_HA_siNC_mnase32_x_chip_rep1.mm9.bw",
                                   TU.DE.mm9.gr, fast = T) %>% unlist() %>% unname()
  
  # Histone marks on TSSs --------------------------------------------------------------------------------------------
  tmp <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/Histone/", ".bw$", full.names = T), 
                  TU.TSS.mm9.gr, fast = T) %>% '*'(width(TU.TSS.mm9.gr) / 1e3)
  tmp <- tmp[, -19]
  colnames(tmp) <- paste0("TSS_", gsub("_ChIP.*|_rep.*|_wt.*|_WT", "\\1", gsub(".*20.*_(H.*).mm9.bw", "\\1", colnames(tmp))))
  TU_feature_table <- cbind(TU_feature_table, tmp)
  
  saveRDS(TU_feature_table, "data/TU_feature_table.RData")
}

# Cluster features --------------------------------------------------------------------------------------------
# feature_data <- log(TU_feature_table[, -c(1:12)] / (rowMedians(as.matrix(TU_feature_table[, -c(1:12)]))))

keep_features <- c("TSS_DNAme.", "GB_DNAme.", "DHS", 
                   "Aff4", "E2f1", "TBP", "TFIID", "Sp1",
                   "cMyc", "Esrrb", "Klf4", "Nanog", "Oct4", "Sox2",
                   "CTCF", "p300", "Hdac1", "Med1", "Brd4", "Yy1", 
                   "Chd1", "Chd2", "Chd9", "Ezh2", "Ring1b", 
                   "HP1a", "LaminB", 
                   "H2AZ", "H2AX", "H33_YFP", 
                   "H2BK20ac", "H3K9ac", "H3K27ac", "H3K56ac", "H3K64ac", "H3K122ac", 
                   "H3K4me1", "TSS_H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K79me2",   
                   "H2AUb", "H2AK119ub1")

feature_types <- c(rep("DNA", 5),
                   rep("General TF", 5), 
                   rep("Pluripotent TF", 6), 
                   rep("Enhancer Activity", 6),
                   rep("Remodeler", 5),
                   rep("Domain", 2),
                   rep("Histone", 3),
                   rep("Acetylation", 6),
                   rep("Methylation", 6),
                   rep("Ubiq.", 2))

# save feature details
if (F) {
  feature_table <- cbind(features = keep_features,
                         feature_types = feature_types[-(1:2)],
                         gene_part = rep("TSS", length(keep_features)),
                         GEO_catalogue = "NA")
  feature_table[feature_table[, "features"] %in% 
                  c("GB_DNAme.", "Aff4", "Yy1", 
                    "Chd1", "Chd2", "Chd9", "Ezh2", "Ring1b", 
                    "HP1a", "LaminB",
                    "H33_YFP", "H2AZ", "H2AX",
                    "H2BK20ac", "H3K9ac", "H3K27ac", "H3K56ac", "H3K64ac", "H3K122ac", 
                    "H3K4me1", "H3K9me3", "H3K27me3", "H3K36me2", "H3K36me3", "H3K79me2",       
                    "H2AUb", "H2AK119ub1"
                  ), "gene_part"] <- "Gene_body"
  write.table(feature_table, "../fig3/data/feature_detail.txt", quote = F, row.names = F, sep = "\t")
}

# prepare features matrix 
feature_data <- cbind(log2(TU_feature_table[, c("CpG_num", "TATA_num")]), 
                      TU_feature_table[, keep_features[1:2]], # not transform for DNA methylation density
                      log2(TU_feature_table[, keep_features[-(1:2)]]) ) %>%
                as.matrix() %>% trim_quantile(0.97)
rownames(feature_data) <- TU_feature_table$gene_id
colnames(feature_data) <- c("CpG_num", "TATA_num", keep_features)

# remove incompleted values and keep TU types
rm.idx <- apply(feature_data, 1, function(x) any(is.na(x) | is.infinite(x))) %>% which()
mRNA.idx <- which(TU.DE.mm9.gr$location == "protein_coding")
mRNA.idx <- mRNA.idx[!mRNA.idx %in% c(rm.idx, non_mRNA.idx)] # 8373
intergenic.idx <- (TU.DE.mm9.gr$location == "intergenic" & 
  (findOverlaps(flank(TU.DE.mm9.gr[TU.DE.mm9.gr$location == "intergenic"],
                      width = 2000, start = F),
                TU.DE.mm9.gr) %>% countSubjectHits() == 0)) %>% which() # remove truncated intergenic TU
intergenic.idx <- intergenic.idx[!intergenic.idx %in% rm.idx] # 7852

get_Rs <- function(response, feature_data, idx, is.cor = F) {
  multi_variance_explained(X = feature_data[idx, ], Y = response[idx], is.cor = is.cor)
}

if (T) {
  # mRNA R2 
  tmp <- apply(log(TU_feature_table[, c("Copy", "Half_life", "TX_RPK")]) %>% trim_quantile(0.95), 2,
               function(x) get_Rs(response = x,
                                  feature_data = feature_data,
                                  idx = mRNA.idx, 
                                  is.cor = F))
  cat(colSums(tmp), "\n")
  # ΣR2    Copy  Half_life    TX_RPK 
  #   0.4651194 0.05179352 0.7240785 
  rownames(tmp) <- c("CpG", "TATA", keep_features)
  colnames(tmp) <- c("Copy", "Half life", "TX")
  res_data <- reshape::melt(tmp)
  colnames(res_data)[3] <- "R2"
  # mRNA correlation
  tmp <- apply(log(TU_feature_table[, c("Copy", "Half_life", "TX_RPK")]) %>% trim_quantile(0.95), 2,
               function(x) get_Rs(response = x, 
                                  feature_data = feature_data,
                                  idx = mRNA.idx,
                                  is.cor = T))
  rownames(tmp) <- c("CpG", "TATA", keep_features)
  colnames(tmp) <- c("Copy", "Half life", "TX")
  res_data <- cbind(res_data, "Correlation" = reshape::melt(tmp)[, 3],
                    TU = "mRNA", Type = feature_types)
  res_data$X1 <- factor(res_data$X1, levels = c("CpG", "TATA", keep_features))
}

if (T) {
  # intergenic RNA R2
  tmp <- apply(log(TU_feature_table[, c("Copy", "Half_life", "TX_RPK")]) %>% trim_quantile(0.95), 2,
               function(x) get_Rs(response = x, 
                                  feature_data = feature_data, 
                                  idx = intergenic.idx, 
                                  is.cor = F))
  cat(colSums(tmp), "\n")
  # ΣR2  Copy  Half_life   TX_RPK 
  # 0.3083363 0.1610435 0.1690043
  rownames(tmp) <- c("CpG", "TATA", keep_features)
  colnames(tmp) <- c("Copy", "Half life", "TX")
  res_data2 <- cbind(reshape::melt(tmp))
  colnames(res_data2)[3] <- "R2"
  # intergenic RNA correlation
  tmp <- apply(log(TU_feature_table[, c("Copy", "Half_life", "TX_RPK")]) %>% trim_quantile(0.95), 2,
               function(x) get_Rs(response = x,
                                  feature_data = feature_data,
                                  idx = intergenic.idx, 
                                  is.cor = T))
  rownames(tmp) <- c("CpG", "TATA", keep_features)
  colnames(tmp) <- c("Copy", "Half life", "TX")
  res_data2 <- cbind(res_data2, "Correlation" = reshape::melt(tmp)[, 3], 
                     TU = "intergenic", Type = feature_types)
}

res_data <- rbind(res_data, res_data2)
res_data$Type <- factor(res_data$Type, levels = unique(feature_types))

ggplot(res_data, aes(x = X1, y = X2)) +
  geom_tile(aes(fill = R2)) + 
  geom_point(aes(colour = Correlation, size = abs(Correlation))) +
  facet_grid(TU~Type, scales = "free", space = "free", switch = "y") +
  xlab("") + ylab("") +
  scale_fill_viridis_c(option = "E", direction = -1, alpha = 0.65) +
  scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkgreen") +
  scale_size(range = c(0.01, 5), guide = F) +
  theme_setting +
  theme(axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 12, hjust = 1),
        strip.placement = "outside",
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 12)) 

ggsave("figs/Fig3_mRNA_ncRNA_copy_half_life_epigenome_explained3.png",
       width = 14, height = 3.2, device = "png")

# group TU TSS marks with Gaussian mixture clustering and vidualise on umap low-dimentional spaces -------------------------------------
keep.idx <- match(res_data[res_data$TU == "mRNA" &
                             res_data$X2 == "TX" & 
                             res_data$R2 > 0, 
                           "X1"], res_data$X1) # 26 signicant features 
feature_data <- feature_data[, keep.idx] %>% scale() 
tmp.idx <- c(mRNA.idx, intergenic.idx)

require(mclust)
mclust_cls <- Mclust(feature_data[tmp.idx, ], G = 6)

# reduce multi-feature matrix dimensionality for plotting
feature_umap <- umap::umap(feature_data[tmp.idx, ])

# make umap colors
umap_dat <- data.frame("UMAP1" = feature_umap$layout[, 1],
                       "UMAP2" = feature_umap$layout[, 2])
umap_dat <- cbind("Cluster" = mclust_cls$classification,
                  umap_dat, 
                  "Type" = TU.DE.mm9.gr$location[tmp.idx],
                  "Copy" = log10(TU_feature_table$Copy[tmp.idx]),
                  "Half life" = log10(TU_feature_table$Half_life[tmp.idx]) %>% trim_quantile(),
                  "Transcription" = log10(TU_feature_table$TX_RPK[tmp.idx]) %>% trim_quantile()
                  )
umap_dat$Type <- as.character(umap_dat$Type)
umap_dat$Type[umap_dat$Type == "protein_coding"] <- "mRNA"
umap_dat$Density <- 0
umap_dat$Density[umap_dat$Type == "intergenic"] <- with(umap_dat[umap_dat$Type == "intergenic", ],
                                                        get_dens(UMAP1, UMAP2, n.grid = 200))
umap_dat$Density[umap_dat$Type != "intergenic"] <- with(umap_dat[umap_dat$Type != "intergenic", ],
                                                        get_dens(UMAP1, UMAP2, n.grid = 200))
# reorder cluster by tx level
umap_dat$Cluster <- factor(umap_dat$Cluster, 
                           levels = umap_dat %>%
                             dplyr::group_by(Cluster) %>%
                             dplyr::summarise(means = mean(Transcription[!is.infinite(Transcription)&!is.na(Transcription)])) %>% 
                             "$"("means") %>% order(decreasing = T) %>% unname())
umap_dat$Cluster <- as.factor(as.numeric(umap_dat$Cluster))

g1 <- umap_dat %>%
  dplyr::mutate("Density" = (Density / max(Density))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = Type, color = Type, alpha = 0.5)) +
  geom_point(pch = 19, size = 1, color = "grey50", stroke = 0.5) +
  geom_point(pch = 19, size = 1, color = "white") +
  geom_point(pch = 19, size = 0.3) +
  scale_color_manual(values = c("tan2", "steelblue")) +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         alpha = FALSE)

for (i in levels(umap_dat$Cluster)) {
  umap_dat$Density[umap_dat$Cluster == i] <- with(umap_dat[umap_dat$Cluster == i, ],
                                                  get_dens(UMAP1, UMAP2, n.grid = 200))
}

g2 <- umap_dat %>%
  # dplyr::mutate("Density" = get_dens(UMAP1, UMAP2, n.grid = 200)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Cluster, alpha = 0.5)) +
  geom_point(pch = 19, size = 1, color = "grey50", stroke = 0.5) +
  geom_point(pch = 19, size = 1, color = "white") +
  geom_point(pch = 19, size = 0.3) +
  scale_color_brewer(type='qual', palette = 6, direction = 1) +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)),
         alpha = FALSE)

g3 <- umap_dat %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate("Density" = get_dens(UMAP1, UMAP2, n.grid = 200)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = `Half life`, alpha = 0.5)) +
  geom_point(pch = 19, size = 1, color = "grey50", stroke = 0.5) +
  geom_point(pch = 19, size = 1, color = "white") +
  geom_point(pch = 19, size = 0.3) +
  scale_color_gradientn(colours = brewer.pal(8, "Blues"), 
                        limits = c(0, 2.5),
                        labels = 10^(0:2),
                        breaks = 0:2) +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0)) +
  guides(alpha = FALSE)

g4 <- umap_dat %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate("Density" = get_dens(UMAP1, UMAP2, n.grid = 200)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Transcription, alpha = 0.5)) +
  geom_point(pch = 19, size = 1, color = "grey50", stroke = 0.5) +
  geom_point(pch = 19, size = 1, color = "white") +
  geom_point(pch = 19, size = 0.3) +
  scale_color_gradientn(name = "Tx RPK",
                        colours = brewer.pal(8, "Greens")) +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0)) +
  guides(alpha = FALSE)

ggsave(plot = cowplot::plot_grid(g1, g2, g3, g4, ncol=2, align = "v"), 
       filename = "../fig3/figs/Fig3_umap_mRNA_intergenic_features3.png",
       width = 7, height = 5, device = "png")

# ---------------------------------------------------------------------------
# show PTM details on the heatmap
PTMs <- c("TSS_H3K4me3", "H2AZ", "H33", "H3K9me3", "H3K27me3", "H2AUb", "H3K27ac", "H3K79me2")
heatmap_dat <- data.frame("Tx RPK" = umap_dat$Transcription, #TU_feature_table$TX_RPK[tmp.idx] %>% log10(),
                          "Half-life" = umap_dat$`Half life`, #TU_feature_table$Half_life[tmp.idx] %>% log10(),
                          "CpG" = log1p(TU_feature_table$CpG_num[tmp.idx]),
                          TU_feature_table[tmp.idx, PTMs] %>% log1p() %>% trim_quantile(0.975))
heatmap_dat <- as.data.frame(scale(heatmap_dat))
rownames(heatmap_dat) <- seq_len(nrow(heatmap_dat))
colnames(heatmap_dat)[-(1:3)] <- gsub("TSS_", "", PTMs)

anno_row <- data.frame("Type" = as.character(umap_dat$Type), 
                       "Cluster" = as.character(umap_dat$Cluster),
                       "ChromHMM" = TU.TSS.mm9.gr$chromHMM[tmp.idx],
                       "Direction" = TU.TSS.mm9.gr$Direction[tmp.idx])
rownames(anno_row) <- rownames(heatmap_dat)
row_ord <- NULL
for (i in seq_along(unique(umap_dat$Cluster))) {
  row_ord <- c(row_ord, 
               rownames(heatmap_dat[umap_dat$Cluster == i, ])[order(anno_row$ChromHMM[anno_row$Cluster == i])])
}
heatmap_dat <- heatmap_dat[row_ord, ]
anno_row <- anno_row[row_ord, ]

annotation_colors <- list(Type = c("steelblue", "tan2"),
                          ChromHMM = colors_20[c(1:3,5,7:9,11:14)],
                          Direction = colors_9[1:2])
names(annotation_colors[[1]]) <- c("mRNA", "intergenic")
names(annotation_colors[[2]]) <- levels(anno_row$ChromHMM)
names(annotation_colors[[3]]) <- c("Bidirection", "Unidirection")

dev.off()
pdf("../fig3/figs/Fig3_mRNA_Intergenic_multi_feature_clustering3.pdf", width = 4.5, height = 4.5)
out <- pheatmap::pheatmap(heatmap_dat, 
                          show_rownames = F,
                          cluster_rows = F,
                          cluster_cols = F,
                          angle_col = 45,
                          fontsize = 8,
                          border_color= T,
                          annotation_row = anno_row[, -2], 
                          na_col = "white",
                          annotation_colors = annotation_colors,
                          gaps_row = cumsum(table(anno_row$Cluster)), 
                          gaps_col = c(2, 4, 6, 9)
                          # main = paste0("mRNA (", length(mRNA.idx),
                          #               "), intergenic RNA (",length(intergenic.idx), ")")
                          )
dev.off()

# specify epiminoc marks levels of each cluster
g_list <- list()
for (i in c("Tx.RPK", "Half.life", "H3K4me3", "H2AZ", "H33", "H3K9me3", 
            "H3K27me3", "H2AUb", "H3K27ac", "H3K79me2")) {
  g_list <- c(g_list, 
              list(data.frame(x = paste0("C", anno_row$Cluster), y = heatmap_dat[, i]) %>% 
                     ggplot(aes(x = x, y = y, fill = x)) +
                     geom_boxplot( notch = T, outlier.size = 0) +
                     xlab("") + ylab("Z-scores") + ggtitle(i) +
                     theme_setting +
                     theme(legend.position = "none")
                   ))
}
ggsave(plot = do.call(grid.arrange, c(g_list, nrow = 2)),
       filename = paste0("Fig3_Boxplot_TU_heatmap_features.png"), 
       path = "../fig3/figs/",
       device = "png", width = 16, height = 6)

# next is "Fig3_tx_changes_vs_histone_changes.R"
