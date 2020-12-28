# Exon Intron TPM
library(org.Mm.eg.db)

# exon / intron ratios ------------------------------------------------------------------------------
gene_ids <- read.table("../data/active_es_genes.txt")[, 1]
# tx_ids <- read.table("../data/active_es_txs.txt")[, 1]
# Exons
exon.parts <- GenomicFeatures::exonicParts(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                           linked.to.single.gene.only = T)
exon.parts <- exon.parts[exon.parts$gene_id %in% gene_ids]
exon.parts <- exon.parts[seqnames(exon.parts) %in% c(1:19, "X", "Y")]
seqlevelsStyle(exon.parts) <- "UCSC"
# Introns
intron.parts <- GenomicFeatures::intronicParts(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                               linked.to.single.gene.only = T)
intron.parts <- intron.parts[intron.parts$gene_id %in% gene_ids]
intron.parts <- intron.parts[seqnames(intron.parts) %in% c(1:19, "X", "Y")]
seqlevelsStyle(intron.parts) <- "UCSC"


# Count bam reads
bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                        pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("SL2i|rep3", bam_files)]

exon_counts <- Rsubread::featureCounts(bam_files,
                                       annot.ext = data.frame(GeneID = exon.parts$gene_id,
                                                              Chr = as.character(seqnames(exon.parts)),
                                                              Start = start(exon.parts),
                                                              End = end(exon.parts),
                                                              Strand = strand(exon.parts)),
                                       isPairedEnd=TRUE, strandSpecific = 1, nthreads = 16)$counts
exon_widths <- aggregate(width(exon.parts)/1e3, list(exon.parts$gene_id), sum)
exon_rpk <- exon_counts / exon_widths[match(rownames(exon_counts), exon_widths$Group.1), 2]

intron_counts <- Rsubread::featureCounts(bam_files, strandSpecific = 1, nonSplitOnly = T,
                                         annot.ext = data.frame(GeneID = intron.parts$gene_id,
                                                                Chr = as.character(seqnames(intron.parts)),
                                                                Start = start(intron.parts),
                                                                End = end(intron.parts),
                                                                Strand = strand(intron.parts)),
                                         isPairedEnd=TRUE, nthreads = 16)$counts
intron_widths <- aggregate(width(intron.parts)/1e3, list(intron.parts$gene_id), sum)
intron_rpk <- intron_counts / intron_widths[match(rownames(intron_counts), intron_widths$Group.1), 2]

gene_ov <- intersect.Vector(rownames(exon_rpk), rownames(intron_rpk))
exon_rpk <- exon_rpk[gene_ov, ]
intron_rpk <- intron_rpk[gene_ov, ]

intron_content_class <- intron_widths$x[match(gene_ov, intron_widths$Group.1)] / 
  exon_widths$x[match(gene_ov, exon_widths$Group.1)]
intron_content_class <- cut(intron_content_class, quantile(intron_content_class, seq(0, 1, 0.2)))
names(intron_content_class) <- gene_ov

intron_length_class <- cut(intron_widths$x[match(gene_ov, intron_widths$Group.1)], 
                           quantile(intron_widths$x, seq(0, 1, 0.2)))
names(intron_length_class) <- gene_ov

# combine replicates after normalization
LRNA.sizefactor <- readRDS("../data/LRNA.sizefactor.RC.RData")
FRNA.sizefactor <- readRDS("../data/FRNA.sizefactor.RC.RData")
sizefactors <- c(`names<-`(FRNA.sizefactor, paste0("FRNA_", names(FRNA.sizefactor))), 
                 `names<-`(LRNA.sizefactor, paste0("LRNA_", names(LRNA.sizefactor))))

colnames(exon_rpk) <- gsub("\\.", "_", gsub(".Aligned.*", "", colnames(exon_rpk)))
colnames(intron_rpk) <- gsub("\\.", "_", gsub(".Aligned.*", "", colnames(intron_rpk)))

exon_rpk <- sweep(exon_rpk, 2, sizefactors[colnames(exon_rpk)], "/") 
exon_rpk <- `colnames<-`(t(apply(exon_rpk, 1, function(x) 
  aggregate(x, list(gsub("_rep.*", "", colnames(exon_rpk))), "sum")$x )),
  unique(gsub("_rep.*", "", colnames(exon_rpk))))

intron_rpk <- sweep(intron_rpk, 2, sizefactors[colnames(intron_rpk)], "/") 
intron_rpk <- `colnames<-`(t(apply(intron_rpk, 1, function(x) 
  aggregate(x, list(gsub("_rep.*", "", colnames(intron_rpk))), "sum")$x )),
  unique(gsub("_rep.*", "", colnames(intron_rpk))))

rm_idx <- apply(cbind(exon_rpk, intron_rpk), 1, function(x) any(x == 0))

intron_exon_ratio <- intron_rpk / exon_rpk # intron spliced rate
intron_exon_ratio <- intron_exon_ratio[!rm_idx & rowMeans(intron_rpk) > 5 & rowMeans(exon_rpk) > 5, ] # keep high expressed genes
# intron_exon_ratio[intron_exon_ratio < 0] <- 0 # remove PCR amplification bias on introns, let intron with larger RPK than exon be unspliced

pre.frac <- intron_exon_ratio[, grep("FRNA", colnames(intron_exon_ratio))]
lab.frac <- intron_exon_ratio[, grep("LRNA", colnames(intron_exon_ratio))]

splicing.rate <- lab.frac / pre.frac

splicing_change_2i <- log(lab.frac[, 5] / lab.frac[, 1])
splicing_change_2i <- splicing_change_2i[splicing_change_2i > (-0.3) & splicing_change_2i < 0.3]

splicing_change_mTORi <- log(lab.frac[, 5] / lab.frac[, 3])
splicing_change_mTORi <- splicing_change_mTORi[splicing_change_mTORi > (-1) & splicing_change_mTORi < 1]


speed_change_2i <- with(gene_body_production_mat[gene_ov, ],
                        (Tx_RPK_SL - Pol2_SL) - (Tx_RPK_2i - Pol2_2i))
names(speed_change_2i) <- gene_ov
# speed_change_2i <- speed_change_2i[speed_change_2i > (- 0.5) & speed_change_2i < 1]

single_variance_explained(splicing_change_2i, 
                          speed_change_2i[match(names(splicing_change_2i), 
                                                names(speed_change_2i))],
                          is.cor = T)

single_variance_explained(splicing_change_mTORi, 
                          speed_change_mTORi[match(names(splicing_change_mTORi),
                                                   names(speed_change_mTORi))], T)


smoothScatter(splicing_change_2i, 
              speed_change_2i[match(names(splicing_change_2i), names(speed_change_2i))], 
              xlim = c(-1, 1))

smoothScatter(splicing_change_mTORi, 
              speed_change_mTORi[match(names(splicing_change_mTORi),
                                       names(speed_change_mTORi))] )


saveRDS(list("pre.frac" = pre.frac, "lab.frac" = lab.frac), 
        "data/splicing_ratio.RData")

# get GRO-seq or PRO-seq coverage on gene exons -----------------------------------------------------
