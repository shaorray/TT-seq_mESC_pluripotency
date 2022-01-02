# Rui Shao
# Jun 2020

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

L_sf = readRDS("../data/LRNA.sizefactor.RData")
F_sf = readRDS("../data/FRNA.sizefactor.RData")

Gene_input = importRanges("/mnt/0E471D453D8EE463/genomeDir/Mus_musculus.GRCm38.79.gtf")
TU_input = importRanges("../data/TU_anno/mm10/TU_filter+LRNA_SL.gff3")

# bam.file.list:
#  - bam.files of different samples (output will average replicates after normalization)
#     - bam.file coverage of a signle replicate (extract and normalise with size factors) 

bam.file.list = list("LRNA_SL" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                           pattern = "LRNA.*(SL_).*bam$", full.names = T),
                     "LRNA_2i" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                            pattern = "LRNA.*(_2i_2d).*bam$", full.names = T),
                     "LRNA_mTORi" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                            pattern = "LRNA.*(mTORi_1d).*bam$", full.names = T),
                     "FRNA_SL" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                            pattern = "FRNA.*(SL_).*bam$", full.names = T),
                     "FRNA_2i" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                            pattern = "FRNA.*(_2i_2d).*bam$", full.names = T),
                     "FRNA_mTORi" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                                            pattern = "FRNA.*(mTORi_1d).*bam$", full.names = T))

size.factor.list = list("LRNA_SL" = L_sf[grep("SL", names(L_sf))],
                        "LRNA_2i" = L_sf[grep("2i_2", names(L_sf))],
                        "LRNA_mTORi" = L_sf[grep("mTORi_1", names(L_sf))],
                        "FRNA_SL" = F_sf[grep("SL", names(F_sf))],
                        "FRNA_2i" = F_sf[grep("2i_2", names(F_sf))],
                        "FRNA_mTORi" = F_sf[grep("mTORi_1", names(F_sf))] )

annotation_list <- vector(mode = "list", length = 2)
names(annotation_list) <- c("Gene", "TU")

annotation_list = list("Gene" = Gene_input[!grepl("Gm", Gene_input$gene_name)],
                       "TU" = TU_input)

cov_color_list = list("LRNA_SL" = c("plus"=colors_20[13], "minus"=add.alpha(colors_20[13],.4)),
                      "LRNA_2i" = c("plus"=colors_20[2], "minus"=add.alpha(colors_20[2],.4)),
                      "LRNA_mTORi" = c("plus"=colors_20[20], "minus"=add.alpha(colors_20[20],.4)),
                      "FRNA_SL" = c("plus"=colors_20[13], "minus"=add.alpha(colors_20[13],.4)),
                      "FRNA_2i" = c("plus"=colors_20[2], "minus"=add.alpha(colors_20[2],.4)),
                      "FRNA_mTORi" = c("plus"=colors_20[20], "minus"=add.alpha(colors_20[20],.4)) )

# Kit
pdf("figs/Fig1_genome_view_kit2.pdf", width = 15, height = 8)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr5', 
                                IRanges(start=73753000, end=75830000), strand = "+"),
             stranded = T, 
             size.factor.list = size.factor.list,
             low_cut = 0.5, 
             log_scale = F, 
             bin_width = 500,
             smoothen = F,
             anno_text_cex = 1,
             anno_text_offset = 1.8,
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)
dev.off()

# Tbx3
pdf("figs/Fig1_genome_view_Tbx3.pdf", width = 6, height = 7)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr5', 
                                IRanges(start = 119650000, end = 119693000), strand = "+"),
             stranded = T, 
             size.factor.list = size.factor.list,
             low_cut = 2.5,
             log_scale = T,
             bin_width = 200,
             smoothen = F,
             anno_text_cex = 1,
             anno_text_offset = 1.8,
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)
dev.off()

# Myc (cMyc) 40 kb window
pdf("figs/Fig1_genome_view_myc.pdf", width = 6, height = 7)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr15', IRanges(start=61977931, end=61977931+40000), strand="+"),
             stranded = T, is.fragment = F,
             size.factor.list = size.factor.list, 
             log_scale = T, 
             smoothen = F, 
             df = 500,
             low_cut = 2.5,
             bin_width = 100, 
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list,
             # anno_box_height = 1.6,
             anno_text_cex = 1,
             anno_text_offset = 2)
dev.off()

# Mycn
pdf("figs/Fig1_genome_view_mycn.pdf", width = 8, height = 10)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr12', IRanges(start=12927193, end=12945119), strand="+"),
             stranded = T, 
             size.factor.list = size.factor.list, 
             log_scale = T, smoothen = F, 
             bin_width = 100, 
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)
dev.off()

# Pou5f1
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr17', IRanges(start=35503263, end=35511121), strand="+"),
             stranded = T, 
             size.factor.list = size.factor.list, 
             log_scale = T, smoothen = T, df = 200,
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)

# ------------------------------------------- Revision ------------------------------------------- #
# (Fig EV7) TTseq / Pol II-S5p, [synthesis] + [velocity] + [pausing] + [termination distance]

bam.file.list = list("LRNA_SL" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                            pattern = "LRNA.*(SL_).*bam$", full.names = T),
                     "Pol2S5p_SL" = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/P1_Pol5p_SL_CTR_ALL.mm9.fltd.bam",
                     
                     "LRNA_2i" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                            pattern = "LRNA.*(_2i_2d).*bam$", full.names = T),
                     "Pol2S5p_2i" = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/P1_Pol5p_2i_CTR_ALL.mm9.fltd.bam",
                     
                     "LRNA_mTORi" = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm9",
                                               pattern = "LRNA.*(mTORi_1d).*bam$", full.names = T),
                     "Pol2S5p_mTORi" = "/mnt/E0767589767560E8/UPPMAX/PUBLISH/bam/P1_Pol5p_SL_INH_ALL.mm9.fltd.bam")

pol2s5p_input_size_factor <- readRDS("../fig3/data/pol2s5p_input_size_factor.RData") %>% "/"(8) # scale up

size.factor.list = list("LRNA_SL" = L_sf[grep("SL_", names(L_sf))],
                        "Pol2S5p_SL" = pol2s5p_input_size_factor[3],
                        
                        "LRNA_2i" = L_sf[grep("2i_2", names(L_sf))],
                        "Pol2S5p_2i" = pol2s5p_input_size_factor[1],
                        
                        "LRNA_mTORi" = L_sf[grep("mTORi_1", names(L_sf))],
                        "Pol2S5p_mTORi" = pol2s5p_input_size_factor[4])


Gene_input = importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
annotation_list = list("Gene" = Gene_input[!grepl("Gm", Gene_input$gene_name) & Gene_input$gene_type == "protein_coding"])

cov_color_list = list("LRNA_SL" = c("plus"=colors_20[13], "minus"=add.alpha(colors_20[13],.4)),
                      "LRNA_2i" = c("plus"=colors_20[2], "minus"=add.alpha(colors_20[2],.4)),
                      "LRNA_mTORi" = c("plus"=colors_20[20], "minus"=add.alpha(colors_20[20],.4)),
                      "Pol2S5p_SL" = c("plus"=add.alpha(colors_20[13],.4)),
                      "Pol2S5p_2i" = c("plus"=add.alpha(colors_20[2],.4)),
                      "Pol2S5p_mTORi" = c("plus"=add.alpha(colors_20[20],.4)) )

# Tcfcp2l1
pdf("../figS6/figs/FigS6_genome_view_Tcfcp2l1.pdf", width = 10, height = 8)
viewCoverage(bam.file.list = bam.file.list,
             interval = GRanges(seqnames = 'chr1', IRanges(start=120510923, end=120605800), strand="+"),
             stranded = F,
             size.factor.list = size.factor.list,
             low_cut = 1,
             log_scale = F,
             bin_width = 300,
             smoothen = F,
             anno_text_cex = 2,
             anno_text_offset = 4,
             cov_color_list = cov_color_list,
             annotation_list = annotation_list)
dev.off()


# Tead4
pdf("../figS6/figs/FigS6_genome_view_Tead4.pdf", width = 10, height = 8)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr6', IRanges(start=128154604, end=128260984), strand="-"),
             stranded = F, 
             size.factor.list = size.factor.list,
             low_cut = 2,
             log_scale = F,
             bin_width = 400,
             smoothen = F,
             anno_text_cex = 2,
             anno_text_offset = 4,
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)
dev.off()


# Etv5
pdf("../figS6/figs/FigS6_genome_view_Etv5.pdf", width = 10, height = 8)
viewCoverage(bam.file.list = bam.file.list, 
             interval = GRanges(seqnames = 'chr16', IRanges(start=22362537, end=22457454), strand="-"),
             stranded = F, 
             size.factor.list = size.factor.list,
             low_cut = 1,
             log_scale = F,
             bin_width = 300,
             smoothen = F,
             anno_text_cex = 2,
             anno_text_offset = 4,
             cov_color_list = cov_color_list, 
             annotation_list = annotation_list)
dev.off()
