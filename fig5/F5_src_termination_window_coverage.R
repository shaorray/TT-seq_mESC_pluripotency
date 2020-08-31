
if (!file.exists("data/last.exon.cov.MINUTE_normed_Pol2S5p.SL.2i.RData"))
{

  cov.last.exon.list = readBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                                                      pattern = "P1.*5p.*ALL.*bam$", full.names = T),
                               intervals = gene.last.exon.gr,
                               pair_end = T,
                               stranded = F,
                               flanks = c(2000, 15000),
                               new_lens = c(50, 50, 200))
  
  input_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                           pattern = "P1.*IN_.*ALL.*bam$", full.names = T)

  sf = SizeFactorCal(.countBam(bam_files = input_files, 
                               intervals = gene.last.exon.gr))
  
  # normalize coverages
  cov.last.exon.list_normed <- list()
  for( i in seq_along(sf) )
    cov.last.exon.list_normed <- c(cov.last.exon.list_normed,
                                  list(cov.last.exon.list[[i]] / sf[i]))
  
  saveRDS(cov.last.exon.list_normed, "data/last.exon.cov.MINUTE_normed_Pol2S5p.SL.2i.RData")
}

if (F) {
  # process coverage on termination windows
  # Pol2S5p
  bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam",
                         pattern = "P1.*5p.*ALL.*bam$", full.names = T)
  
  Pol2S5p.MINUTE.terWin.cov = .coverBam(bam_files, paired.end = "extend",
                                        intervals = terWindow, is.resize = F)
  saveRDS(Pol2S5p.MINUTE.terWin.cov, "data/Pol2S5p.MINUTE.terWin.cov.RData")
  
  # TT-seq
  bam.file.iist = list("LRNA_SL" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                              pattern = "LRNA.*(SL_).*bam$", full.names = T),
                       "LRNA_2i" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                              pattern = "LRNA.*(_2i_2d).*bam$", full.names = T),
                       "LRNA_mTORi" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                                 pattern = "LRNA.*(mTORi_2d).*bam$", full.names = T))
  
  TTseq.terWin.cov = lapply(bam.file.iist, .coverBam, 
                            intervals = terWindow, paired.end = "extend", is.resize = F)
  
  L_sf = readRDS("../data/LRNA.sizefactor.RData")
  size.factor.list = list("LRNA_SL" = L_sf[grep("SL_", names(L_sf))],
                          "LRNA_2i" = L_sf[grep("^2i_2", names(L_sf))],
                          "LRNA_mTORi" = L_sf[grep("mTORi_2", names(L_sf))])
  # normalize with size factors, and combine replicates
  TTseq.terWin.cov_norm <- vector(mode = "list", length = length(TTseq.terWin.cov))
  names(TTseq.terWin.cov_norm) <- names(TTseq.terWin.cov)
  for (sample in names(TTseq.terWin.cov)) { 
    TTseq.terWin.cov_norm[[sample]] <- 
      sapply(seq_along(TTseq.terWin.cov[[sample]][[1]]), 
             function(i) 
               rowMeans(sapply(seq_along(TTseq.terWin.cov[[sample]]), function(x) TTseq.terWin.cov[[sample]][[x]][i] / size.factor.list[[sample]][x]))
      )
  }
  saveRDS(TTseq.terWin.cov_norm, "data/TTseq.terWin.cov_norm.RData")
  
  # other nascent RNA-seq
  bam.file.iist = list("Jan_TTseq" = list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Jan_TTseq_TX1072_1/bam/",
                                                     pattern = "out.bam$", full.names = T),
                       "Jesse_PROseq" = c(list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Jesse_PROSeq_ControlClone2_Rep1/bam",
                                                   pattern = "out.bam$", full.names = T),
                                          list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Jesse_PROSeq_ControlClone2_Rep2/bam",
                                                     pattern = "out.bam$", full.names = T)),
                       "Marta_PROseq" = c(list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Marta_ESC_1_PRO-Seq_ctrl/bam",
                                                     pattern = "out.bam$", full.names = T),
                                          list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Marta_ESC_2_PRO-Seq_ctrl/bam",
                                                     pattern = "out.bam$", full.names = T))
                       # "Williams_C2_2i_GROseq" = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2015_Williams_Bl6_mESC_C2cells_GROseq",
                       #                                      pattern = "2i_GRO.*bam$", full.names = T),
                       # "Flynn_GROseq" = c(list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Flynn_GRO-seq-mES-12hrKD_Control1/bam",
                       #                               pattern = "out.bam$", full.names = T),
                       #                    list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Flynn_GRO-seq-mES-12hrKD_Control2/bam",
                       #                               pattern = "out.bam$", full.names = T),
                       #                    list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Flynn_GRO-seq-mES-12hrKD_Control3/bam",
                       #                               pattern = "out.bam$", full.names = T),
                       #                    list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2016_Flynn_GRO-seq-mES-12hrKD_Control4/bam",
                       #                               pattern = "out.bam$", full.names = T)),
                       # "Dorighi_GROseq" = list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2017_Dorighi_Gro-seq_WT/bam",
                       #                               pattern = "out.bam$", full.names = T),
                       # "Mylonas_NETseq" = c(list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Mylonas_NET-seq_Control_rep1/bam",
                       #                                 pattern = "out.bam$", full.names = T),
                       #                      list.files("/mnt/0E471D453D8EE463/nascent_RNA_mm9/2018_Mylonas_NET-seq_Control_rep2/bam",
                       #                                 pattern = "out.bam$", full.names = T))
                       )
  
  nascent.terWin.cov = lapply(bam.file.iist, .coverBam, 
                              intervals = terWindow, is.resize = F)
  nascent.terWin.cov_norm <- vector(mode = "list", length = length(nascent.terWin.cov))
  names(nascent.terWin.cov_norm) <- names(nascent.terWin.cov)
  for (sample in names(nascent.terWin.cov)) {
    nascent.terWin.cov_norm[[sample]] <- 
      sapply(seq_along(nascent.terWin.cov[[sample]][[1]]), 
             function(i) 
               rowMeans(sapply(seq_along(nascent.terWin.cov[[sample]]), function(x) nascent.terWin.cov[[sample]][[x]][i]))
      )
  }
  saveRDS(nascent.terWin.cov_norm, "data/nascent.terWin.cov_norm.RData")
  
  # other Pol II ChIP
  bam.file.iist = list("Ferrai_Pol_II_S5p" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2017_Ferrai_ChIP_RNAPIIS5p_ESC/bam",
                                                        pattern = "Ferrai.*bam$", full.names = T),
                       "Ferrai_Pol_II_S7p" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2017_Ferrai_ChIP_RNAPIIS7p_ESC/bam",
                                                        pattern = "Ferrai.*bam$", full.names = T),
                       "Flynn_Pol_II_DMSO" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2016_Flynn_ChIPseq_mES_Pol_II_DMSO_JQ1/",
                                                        pattern = "DMSO.*bam$", full.names = T),
                       "Flynn_Pol_II_JQ1" = list.files("/mnt/0E471D453D8EE463/GEO_pol2/2016_Flynn_ChIPseq_mES_Pol_II_DMSO_JQ1/",
                                                       pattern = "JQ1.*bam$", full.names = T))
  pol2.terWin.cov = lapply(bam.file.iist, .coverBam, 
                           intervals = terWindow, is.resize = F)
  pol2.terWin.cov_norm <- vector(mode = "list", length = length(pol2.terWin.cov))
  names(pol2.terWin.cov_norm) <- names(pol2.terWin.cov)
  for (sample in names(pol2.terWin.cov)) {
    pol2.terWin.cov_norm[[sample]] <- 
      sapply(seq_along(pol2.terWin.cov[[sample]][[1]]), 
             function(i) 
               rowMeans(sapply(seq_along(pol2.terWin.cov[[sample]]), function(x) pol2.terWin.cov[[sample]][[x]][i]))
      )
  }
  saveRDS(pol2.terWin.cov_norm, "data/pol2.terWin.cov_norm.RData")
}

