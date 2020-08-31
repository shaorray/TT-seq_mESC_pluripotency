# Exon Intron TPM
library(dplyr)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)

gene_ids <- read.table("../data/active_es_genes.txt")[, 1]
tx_ids <- read.table("../data/active_es_txs.txt")[, 1]
# Exons
exon.parts <- GenomicFeatures::exonicParts(EnsDb.Mmusculus.v79,
                                           linked.to.single.gene.only = T)
exon.parts <- exon.parts[sum(exon.parts$gene_id %in% gene_ids) > 0 &
                           sum(exon.parts$tx_id %in% tx_ids) > 0]

# Introns
intron.parts <- GenomicFeatures::intronicParts(EnsDb.Mmusculus.v79,
                                               linked.to.single.gene.only = T)
intron.parts <- intron.parts[sum(intron.parts$gene_id %in% gene_ids) > 0 &
                           sum(intron.parts$tx_id %in% tx_ids) > 0]

# Count bam reads
bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10",
                        pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("SL2i|rep3", bam_files)]
exon_counts <- Rsubread::featureCounts(bam_files,
                                       annot.ext = data.frame(GeneID = exon.parts$gene_id,
                                                              Chr = paste0("chr", seqnames(exon.parts)),
                                                              Start = start(exon.parts),
                                                              End = end(exon.parts),
                                                              Strand = strand(exon.parts)),
                                       isPairedEnd=TRUE, strandSpecific = 1, nthreads = 8)
exon_counts <- exon_counts$counts
exon_counts <- exon_counts[order(rownames(exon_counts)), ]
exon_rpk <- exon_counts /
  aggregate(width(exon.parts)/1e3, list(exon.parts$gene_id), sum)[, 2]


intron_counts <- Rsubread::featureCounts(bam_files, strandSpecific = 1, nonSplitOnly = T,
                                         annot.ext = data.frame(GeneID = intron.parts$gene_id,
                                                                Chr = paste0("chr", seqnames(intron.parts)),
                                                                Start = start(intron.parts),
                                                                End = end(intron.parts),
                                                                Strand = strand(intron.parts)),
                                         isPairedEnd=TRUE, nthreads = 8)$counts

intron_counts <- intron_counts[order(rownames(intron_counts)), ]
intron_rpk <- intron_counts /
  aggregate(width(intron.parts)/1e3, list(intron.parts$gene_id), sum)[, 2]

gene_ov <- intersect.Vector(rownames(exon_rpk), rownames(intron_rpk))
intron_exon_ratio <- 1 - intron_rpk[gene_ov, ] / exon_rpk[gene_ov, ] # intron spliced proportion
high_expr_gene <- rowMeans(intron_rpk[gene_ov, ]) > 10 & rowMeans(exon_rpk[gene_ov, ]) > 10
intron_exon_ratio <- intron_exon_ratio[high_expr_gene, ]

rm_idx <- apply(intron_exon_ratio, 1, function(x) any(is.na(x) | is.infinite(x) | x == 0))

pre.frac <- intron_exon_ratio[!rm_idx, grep("FRNA", colnames(intron_exon_ratio))] %>% log10
lab.frac <- intron_exon_ratio[!rm_idx, grep("LRNA", colnames(intron_exon_ratio))] %>% log10

saveRDS(list("pre.frac"=pre.frac, "lab.frac"=lab.frac), 
        "data/splicing_ratio.RData")

