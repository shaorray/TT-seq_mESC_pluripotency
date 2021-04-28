fetch_TU_feature <- function(TU.gr, is.fast = T) {
  # get read density from features that have positive R2 in multivariate explanation of mRNA TX
  TSS.gr <- promoters(TU.gr, upstream = 500, downstream = 500)
  TU_feature_table <- data.frame(idx = seq_along(TU.gr))
  # DNA elements on TSS --------------------------------------------------------------------------------------------
  TU_feature_table$CpG_num <- 
    interval_pattern_hits(intervals = promoters(TU.gr, 
                                                upstream = 1000,
                                                downstream = 50),
                          pattern = "CG", 
                          which_genome = "mm9")
  
  TU_feature_table$TATA_num <-
    interval_pattern_hits(intervals = promoters(TU.gr, 
                                                upstream = 1000,
                                                downstream = 50),
                          pattern = "TATA", # using c("TATAAAA", "TATATAA", "TATAAAT", "TATATAT") gives similar results
                          which_genome = "mm9")
  # DNA methylation -----------------------------------------------------------------------------------------------------
  TU_feature_table$TSS_DNAme. <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/CpG/GSM1127953_DNA_CpG_E14_serum_rep.bw", 
                                          TSS.gr, fast = is.fast) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  TU_feature_table$GB_DNAme. <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/CpG/GSM1127953_DNA_CpG_E14_serum_rep.bw", 
                                         TU.gr, fast = is.fast) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  # Chromatin accessibility --------------------------------------------------------------------------------------------
  TU_feature_table$DHS <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/DHS", full.names = T), 
                                   TSS.gr, fast = is.fast) %>% rowMeans() %>% 
    unlist() %>% unname()
  
  TU_feature_table$FAIRE <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/FAIRE", "*", full.names = T), 
                                     TSS.gr, fast = is.fast) %>% rowMeans() %>% 
    unlist() %>% unname()

  # m6A RNA modification --------------------------------------------------------------------------------------------
  TU_feature_table$m6A <- (.countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/m6A", "SySy_ChrMeRIP", full.names = T), 
                                    TU.gr, fast = is.fast) %>% rowMeans() %>% 
                             unlist() %>% unname()) / 
    (.countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/m6A", "SySy_Input", full.names = T), 
              TU.gr, fast = is.fast) %>% rowMeans() %>% 
       unlist() %>% unname())
  
  # single-strand DNA --------------------------------------------------------------------------------------------
  TU_feature_table$KAS <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/ssDNA", "*", full.names = T), 
                                   TSS.gr, fast = is.fast) %>% rowMeans() %>% 
    unlist() %>% unname()

  # General transcription factors / chromatin remodelers on TSS -----------------------------------------------------------------------
  TU_feature_table$Aff4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2011_Lin_AFF4_0hr_ChIPSeq.mm9.bw",
                                    TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$E2f1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2008_Chen_E2f1.mm9.bw",
                                    TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$TBP <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2015_Langer_TBP_ESwt.mm9.bw", 
                                   TSS.gr, fast = is.fast) %>% 
    unlist() %>% unname()
  
  TU_feature_table$TFIID <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2012_Ku_TFIID.mm9.bw",
                                     TSS.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$Sp1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                              "Sp1_rep", full.names = T), 
                                   TSS.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$cMyc <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_cMyc_ChIP_seq.mm9.bw", 
                                    TU.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Esrrb <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Esrrb_ChIP-seq.mm9.bw", 
                                     TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Klf4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Klf4_ChIP_seq.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Nanog <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Nanog_ChIP-seq.mm9.bw", 
                                     TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Oct4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Oct4_ChIP_seq.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Sox2 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Sox2_ChIP_seq.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$CTCF <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2019_Atlasi_GSE92407", 
                                               "CTCF-ChIP-serum", full.names = T), 
                                    TSS.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$p300 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_p300_ChIP-Seq.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Hdac1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2016_Chronis_ESC_Hdac1_ChIP-Seq.mm9.bw", 
                                     TSS.gr, fast = is.fast) %>%
    unlist() %>% unname()
  
  TU_feature_table$Med1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TF/2018_Sabari_MED1_Control.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>% 
    unlist() %>% unname()
  
  TU_feature_table$Brd4 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/TX/2018_Sabari_BRD4_Control.mm9.bw", 
                                    TSS.gr, fast = is.fast) %>% 
    unlist() %>% unname()
  
  TU_feature_table$Yy1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2019_Atlasi_GSE92407", 
                                              "Yy1-ChIP-serum", full.names = T), 
                                   TU.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$Chd1 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                               "Chd1", full.names = T),
                                    TU.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$Chd2 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                               "Chd2", full.names = T),
                                    TU.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$Chd9 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/2016_Hmitou", 
                                               "Chd9", full.names = T),
                                    TU.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$Ezh2 <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/TF", 
                                               "EZH2", full.names = T), 
                                    TSS.gr, fast = is.fast) %>% rowMeans()
  
  TU_feature_table$Ring1b <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/domain/GSM2393579_ES_Ring1b.density.bigWig",
                                      TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$HP1a <- .countBW(list.files("/mnt/0E471D453D8EE463/GEO_bw/domain", "HP1a", full.names = T),
                                    TU.gr, fast = is.fast) %>% rowMeans() %>% unname()
  
  TU_feature_table$LaminB <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/domain/2017_Poleshko_ESC_LaminB_ChIP-seq.mm9.bw",
                                      TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H2AZ <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Surface_H2AZ_WT_ChIP_Seq.mm9.bw",
                                    TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H2AX <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Wu_ChIP_WT_H2AX.mm9.bw",
                                    TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H33_YFP <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/GSM2582412_ESC_H3.3WT_YFP.bigwig",
                                       TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H2BK20ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2015_Kumar_mESC_H2BK20ac_ChIPSeq.mm9.bw",
                                        TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K9ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K9ac_ChIP-Seq.mm9.bw",
                                      TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K27ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K27ac_ChIP-Seq.mm9.bw",
                                       TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K56ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2015_Etchegaray_H3K56ac_ChIP-Seq_Sirt6_WT.mm9.bw",
                                       TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K64ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2019_Martire_ESC_H3K64ac.mm9.bw",
                                       TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K122ac <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2019_Martire_ESC_H3K122ac.mm9.bw",
                                        TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$TSS_H3K4me1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K4me1_ChIP-Seq.mm9.bw",
                                           TSS.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$TSS_H3K4me3 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K4me3_ChIP-Seq.mm9.bw",
                                           TSS.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K9me3 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K9me3_ChIP-Seq.mm9.bw",
                                       TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K27me3 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K27me3_ChIP-Seq.mm9.bw",
                                        TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K36me3 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K36me3_ChIP-Seq.mm9.bw",
                                        TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H3K79me2 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K79me2_ChIP-Seq.mm9.bw",
                                        TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H2AUb <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Kundu_ESC_H2AUb_ChIP.mm9.bw",
                                     TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  TU_feature_table$H2AK119ub1 <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2018_Yao_ES_H2AK119ub1.mm9.bw",
                                          TU.gr, fast = is.fast) %>% unlist() %>% unname()
  
  return(TU_feature_table)
}
