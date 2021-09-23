# Rui Shao 2020 May

# -----------------------------------------------------------------------------------------------------
# This part includes:
#     1. TU transcription level correlation by TU positioning
#     2. intergenic TU occurrence in gene neighborhood, correlation with neighbored genes
#     3. scRNA covariance test of gene neighbor co-expression by positioning
# -----------------------------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# -----------------------------------------------------------------------------------------------------
# load references
TU.DE.mm10.gr <- readRDS("../fig1/data/TU.DE.mm10.RData")
TU.DE.mm10.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$baseMean_LRNA > 0]

gene.gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
gene.gr <- gene.gr[gene.gr$gene_biotype == 'protein_coding']
gene.gr <-  gene.gr[width(gene.gr) < 2500000]

mtch <- findOverlaps(TU.DE.mm10.gr, gene.gr)
TU.DE.mm10.gr$gene_id <- NA
TU.DE.mm10.gr$gene_id[queryHits(mtch)] <- gene.gr$gene_id[subjectHits(mtch)]

# -----------------------------------------------------------------------------------------------------
# gene-gene relative positions
TU.coding.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$location == "protein_coding"]
names(TU.coding.gr) <- TU.coding.gr$gene_id
TU.coding.promoters <- promoters(TU.coding.gr, upstream = 2000, downstream = 0)
levels(strand(TU.coding.promoters)) <- c("-", "+", "*")
TU.coding.ds <- flank(TU.coding.gr, width = 2000, start = F)
levels(strand(TU.coding.ds)) <- c("-", "+", "*")

if (T) {
  # get matches
  #   ds downstream sense
  #   ua upstream antisense
  #   con convergent
  ds_mtch <- findOverlaps(query = TU.coding.gr, 
                          promoters(TU.coding.gr, upstream = 5000, downstream = 0))
  ua_mtch <- findOverlaps(query = TU.coding.gr, subject = TU.coding.promoters)
  con_mtch <- findOverlaps(query = TU.coding.gr, subject = TU.coding.ds)
  
  # downstream sense
  q.hits <- queryHits(ds_mtch)
  s.hits <- subjectHits(ds_mtch)
  
  multi <- ds_mtch[q.hits %in% q.hits[duplicated(q.hits)] | s.hits %in% s.hits[duplicated(s.hits)] ] # 3 or more cis-genes
  multi <- multi[abs(queryHits(multi) - subjectHits(multi)) == 1, ]
  
  uniq <- ds_mtch[!q.hits %in% q.hits[duplicated(q.hits)] & !s.hits %in% s.hits[duplicated(s.hits)] ]
  dse_genes <- cbind("Upstream" = TU.coding.gr$gene_id[c(queryHits(multi), queryHits(uniq))],
                     "Downstream" = TU.coding.gr$gene_id[c(subjectHits(multi), subjectHits(uniq))])
  
  # divergent
  p.embeded <- subjectHits(ua_mtch) %in%
    (findOverlaps(TU.coding.promoters, TU.coding.gr, type = "within") %>% queryHits())
  ua_mtch <- ua_mtch[!p.embeded]
  
  q.hits <- queryHits(ua_mtch)
  s.hits <- subjectHits(ua_mtch)
  multi <- ua_mtch[q.hits %in% q.hits[duplicated(q.hits)] | s.hits %in% s.hits[duplicated(s.hits)] ]
  ms.hits <- subjectHits(multi)
  p.multi <- ms.hits[duplicated(ms.hits)]
  multi.mtch <- multi[ms.hits %in% p.multi]
  multi.strd <- TU.coding.gr[p.multi] %>% strand %>% as.character
  
  p.closest <- sapply(seq_along(p.multi), function(x)
  {
    hits <- split(queryHits(multi.mtch), subjectHits(multi.mtch))[[x]]
    if(multi.strd[x] == "+")
    {
      min(hits)
    } else {
      max(hits)
    }
  }
  )
  
  uniq <- ua_mtch[!q.hits %in% q.hits[duplicated(q.hits)] & !s.hits %in% s.hits[duplicated(s.hits)] ]
  div_genes <- cbind("Promoter" = TU.coding.gr$gene_id[c(p.multi, queryHits(uniq))],
                     "Gene" = TU.coding.gr$gene_id[c(p.closest, subjectHits(uniq))])
  
  # convergent
  q.hits <- queryHits(con_mtch)
  s.hits <- subjectHits(con_mtch)
  multi <- con_mtch[q.hits %in% q.hits[duplicated(q.hits)] ]
  
  f.multi <- unique(queryHits(multi))
  multi.strd <- TU.coding.gr[f.multi] %>% strand %>% as.character
  
  f.closest <- sapply(seq_along(f.multi), function(x)
  {
    hits <- split(subjectHits(multi), queryHits(multi))[[x]]
    if(multi.strd[x] == "+")
    {
      min(hits)
    } else {
      max(hits)
    }
  }
  )
  
  uniq <- con_mtch[ !q.hits %in% q.hits[duplicated(q.hits)] ]
  con_genes <- cbind("Forward" = TU.coding.gr$gene_id[c(f.multi, queryHits(uniq))],
                     "Reverse" = TU.coding.gr$gene_id[c(f.closest, subjectHits(uniq))])
  
  neighbor_gene_list <- list("downstream" = dse_genes, "divergent" = div_genes, "convergent" = con_genes)
  saveRDS(neighbor_gene_list, "data/neighbor_gene_list.RData")
}

# -----------------------------------------------------------------------------------------------------
# gene-intergenic TU neighboring
TU.intergenic.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$location == "intergenic"]
TU.intergenic.gr$id <- seq_along(TU.intergenic.gr)

gene.intergenic.pairs <- foreach (i = seqlevels(TU.coding.gr), .combine = rbind) %dopar% {
  # extract intergenic neighbors for each chromosome
  genes.tmp = TU.coding.gr[seqnames(TU.coding.gr) == i]
  genes.tmp = genes.tmp[order(start(genes.tmp))]
  TU.tmp = TU.intergenic.gr[seqnames(TU.intergenic.gr) == i]
  TU.strands = as.character(strand(TU.tmp)) == "+"
  
  out_table <- NULL
  for (j in seq_along(genes.tmp)) {
    gene.strand = as.character(strand(genes.tmp[j])) == "+"
    gaps = start(TU.tmp) - start(genes.tmp[j])
    
    same.strand = (TU.strands == gene.strand)
    
    if (j == 1) {
      up.tu =  .ifelse(gene.strand, gaps < 0, start(TU.tmp) < start(genes.tmp[2]) & gaps > 0)
      down.tu = .ifelse(!gene.strand, gaps < 0, start(TU.tmp) < start(genes.tmp[2]) & gaps > 0)
    } else if (j == length(genes.tmp)) {
      up.tu =  .ifelse(gene.strand, 
                       start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                       gaps > 0)
      down.tu = .ifelse(!gene.strand,
                        start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                        gaps > 0)
    } else {
      up.tu =  .ifelse(gene.strand, 
                       gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]), 
                       gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
      down.tu = .ifelse(!gene.strand, 
                        gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]),
                        gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
    }
    up.s = same.strand & up.tu
    up.as = !same.strand & up.tu
    down.s = same.strand & down.tu
    down.as = !same.strand & down.tu
    
    out_table <- rbind(out_table,
                       data.frame("gene_id" = rep(genes.tmp[j]$gene_id, sum(up.s + up.as + down.s + down.as)),
                                  "strand" = rep(strand(genes.tmp[j]), sum(up.s + up.as + down.s + down.as)),
                                  "type" = c(rep("up.sense", sum(up.s)), 
                                             rep("up.antisense", sum(up.as)),
                                             rep("down.sense", sum(down.s)),
                                             rep("down.antisense", sum(down.as)) ),
                                  "gap" = c(gaps[up.s], gaps[up.as], gaps[down.s], gaps[down.as]),
                                  "id" = c(TU.tmp$id[up.s], TU.tmp$id[up.as], TU.tmp$id[down.s], TU.tmp$id[down.as])))
  }
  out_table
}

# gap of TU TSS to gene boundary
gaps <- start(promoters(TU.intergenic.gr[gene.intergenic.pairs$id], upstream = 0, downstream = 0)) - 
  cbind(start(TU.coding.gr[gene.intergenic.pairs$gene_id]),
        end(TU.coding.gr[gene.intergenic.pairs$gene_id]))
idx <- ifelse(gene.intergenic.pairs$strand == "+", 
              ifelse(grepl("up", gene.intergenic.pairs$type), 1, 2), 
              ifelse(grepl("up", gene.intergenic.pairs$type), 2, 1))

gene.intergenic.pairs$gap_boundary <- 
  sapply(seq_along(idx),
         function(x) gaps[x, idx[x]]
         ) * ifelse(gene.intergenic.pairs$strand == "+", 1, -1) / 1e3

tmp_table <- cbind(mcols(TU.coding.gr[gene.intergenic.pairs$gene_id])[, 18:22], 
                   mcols(TU.intergenic.gr[gene.intergenic.pairs$id])[, 18:22])
colnames(tmp_table) <- paste0(c(rep("Gene_", 5), rep("TU_", 5)), colnames(tmp_table))
gene.intergenic.pairs <- cbind(gene.intergenic.pairs, tmp_table)
rm(tmp_table)

# show co-expression ------------------------------------------------------------------
gene.intergenic.pairs <-
  gene.intergenic.pairs[!apply(gene.intergenic.pairs[, c(7, 12)], 1, function(x)
    any(is.na(x))), ]

gap_breaks <- c(min(gene.intergenic.pairs$gap_boundary), 
                c(-400, -200, -100, -75, -50, -30, -20, -10, 0, 10, 20, 30, 50, 75, 100, 200, 400),# * 1000,
                max(gene.intergenic.pairs$gap_boundary))

gene.intergenic.pairs$gap_class <- cut(gene.intergenic.pairs$gap_boundary, 
                                       breaks = gap_breaks, 
                                       labels = 1:18) %>% as.numeric()
gene.intergenic.pairs <- gene.intergenic.pairs[!is.na(gene.intergenic.pairs$gap_class), ]

# plot intergenic TU coverage
tmp.intergenic.pairs <- gene.intergenic.pairs[abs(gene.intergenic.pairs$gap_boundary) < 1e2, ]
tmp.gr <- TU.intergenic.gr[tmp.intergenic.pairs$id]

as.idx <- grepl("antisense", tmp.intergenic.pairs$type)
down.idx <- grepl("down", tmp.intergenic.pairs$type)

tmp.cov.gr <- IRanges(start = tmp.intergenic.pairs$gap_boundary + 1e2 - ifelse(down.idx, width(tmp.gr) / 1e3, 0), 
                      width = width(tmp.gr) / 1e3)

s.cov <- coverage(tmp.cov.gr[!as.idx], width = 2e2) %>% spline(x = seq_len(2e2), n = 200) %>% "$"(y)
as.cov <- coverage(tmp.cov.gr[as.idx], width = 2e2) %>% spline(x = seq_len(2e2), n = 200) %>% "$"(y)

insert_NA <- function(x, p = 100, n_0 = 30) c(x[seq_len(p)], rep(NA, n_0), x[seq_len(p) + p])
dat <- data.frame(Pos = rep(1:230, 2),
                  Cov = c(s.cov %>% insert_NA(), 
                          as.cov %>% insert_NA()),
                  Direction = c(rep("Same", 230), rep("Opposite", 230)))
dat$Direction <- factor(dat$Direction, levels = c("Same", "Opposite"))

ggplot(dat, aes(x = Pos, y = Cov)) +
  geom_rect(aes(xmin = 100, xmax = 130, ymin = -Inf, ymax = 0),
            fill = "#0000B2", linetype = 0) +
  geom_rect(aes(xmin = 100, xmax = 130, ymin = 0, ymax = Inf),
            fill = "grey85", linetype = 0) +
  geom_line(aes(color = Direction), lwd = 1) +
  scale_color_manual(name = "Intergenic", values = colors_9[7:6]) +
  scale_x_continuous(name = "\nGene distance (kb)\n",
                     breaks = c(0, 50, 100, 130, 180, 230), 
                     labels = c("-100", "-50", "5'", "3'", "50", "100")) +
  ylab("TU coverage") +
  theme_setting

ggsave(filename = "FigS2_TU_gene_neighbor_coverage_direction.png", 
       path = "../figS2/figs",
       device = "png", width = 6, height = 3.5 )

# enhancer coverage
enhancer_mm10 <- readRDS("../data/enhancer_types.mm10.RData")
enh.idx <- findOverlaps(tmp.gr, enhancer_mm10) %>% countQueryHits() > 0

enh.cov <- coverage(tmp.cov.gr[enh.idx], width = 2e2) %>%
  # smooth.spline(x = seq_len(2e5), df = 20) %>% "$"(y) %>%
  spline(x = seq_len(2e2), n = 200) %>% "$"(y) 
nenh.cov <- coverage(tmp.cov.gr[!enh.idx], width = 2e2) %>% 
  # smooth.spline(x = seq_len(2e5), df = 20) %>% "$"(y) %>%
  spline(x = seq_len(2e2), n = 200) %>% "$"(y)

dat <- data.frame(Pos = rep(1:230, 2),
                  Cov = c(enh.cov %>% insert_NA(), 
                          nenh.cov %>% insert_NA()),
                  Type = c(rep("Enhancer", 230), rep("Other", 230)))
dat$Type <- factor(dat$Type, levels = c("Enhancer", "Other"))

ggplot(dat, aes(x = Pos, y = Cov)) +
  geom_rect(aes(xmin = 100, xmax = 130, ymin = -Inf, ymax = 0),
            fill = "#0000B2", linetype = 0) +
  geom_rect(aes(xmin = 100, xmax = 130, ymin = 0, ymax = Inf),
            fill = "grey85", linetype = 0) +
  geom_line(aes(color = Type), lwd = 1) +
  scale_color_manual(name = "Intergenic", values = colors_9[2:3]) +
  scale_x_continuous(name = "\nGene distance (kb)\n",
                     breaks = c(0, 50, 100, 130, 180, 230), 
                     labels = c("-100", "-50", "5'", "3'", "50", "100")) +
  ylab("TU occurence") +
  theme_setting

ggsave(filename = "FigS2_TU_gene_neighbor_coverage_enhancer.png", 
       path = "../figS2/figs",
       device = "png", width = 6, height = 3.5 )

# plot by TU log2FC correlation
mat <- NULL
for( i in 1:18 ) {
  tmp.pairs <- gene.intergenic.pairs[as.numeric(gene.intergenic.pairs$gap_class) == i, ]
  as.idx <- grepl("antisense", tmp.pairs$type)
  
  enh.idx <- findOverlaps(TU.intergenic.gr[tmp.pairs$id], 
                          subject = enhancer_mm10 ) %>%
             countQueryHits() > 0
  
  tmp.pairs <- tmp.pairs[, grep("log2", colnames(gene.intergenic.pairs))] %>% as.matrix()
  mat <- rbind(mat, c("Pos" = i, 
                      "cor_2i_s" = cor(tmp.pairs[!as.idx, 1], tmp.pairs[!as.idx, 3]),
                      "cor_mTORi_s" = cor(tmp.pairs[!as.idx, 2], tmp.pairs[!as.idx, 4]),
                      "cor_2i_as" = cor(tmp.pairs[as.idx, 1], tmp.pairs[as.idx, 3]),
                      "cor_mTORi_as" = cor(tmp.pairs[as.idx, 2], tmp.pairs[as.idx, 4]),
                      
                      "cor_2i_enh" = cor(tmp.pairs[enh.idx, 1], tmp.pairs[enh.idx, 3]),
                      "cor_mTORi_enh" = cor(tmp.pairs[enh.idx, 2], tmp.pairs[enh.idx, 4]),
                      "cor_2i_nenh" = cor(tmp.pairs[!enh.idx, 1], tmp.pairs[!enh.idx, 3]),
                      "cor_mTORi_nenh" = cor(tmp.pairs[!enh.idx, 2], tmp.pairs[!enh.idx, 4])
                      ))
}
mat[10:18, 1] <- (10:18) + 3 # add gene box
mat <- rbind(mat[1:9, ], matrix(c(10:12, rep(NA, 24)), nrow = 3), mat[10:18, ])


dat_mat <- reshape2::melt(data = as.data.frame(mat[, 1:5]), id = "Pos")
dat_mat$Sample <- rep(c(rep("2i", 21), rep("mTORi", 21)), 2)
dat_mat$Direction <- factor(c(rep("Sense", 21*2), rep("Antisense", 21*2)), 
                        levels = c("Sense", "Antisense"))
dat_mat$line_group <- factor(cumsum(is.na(dat_mat$value)))

g1 <- ggplot(dat_mat, aes(x = Pos, y = value, color = Direction, group = line_group)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey85") +
  # geom_vline(xintercept = 10, col = "grey85") +
  geom_rect(aes(xmin = 10, xmax = 12, ymin = -Inf, ymax = 0),
            fill = "#0000B2", linetype = 0) +
  geom_rect(aes(xmin = 9, xmax = 13, ymin = 0, ymax = Inf),
            fill = "grey95", linetype = 0) +
  geom_rect(aes(xmin = 10, xmax = 12, ymin = 0, ymax = Inf),
            fill = "grey75", linetype = 0) +
  geom_point() + 
  # geom_smooth(method = "gam", span = 0.45, na.rm = FALSE, se = FALSE) +
  geom_line() +
  
  facet_grid(Sample~.) + 
  scale_x_continuous(name = "\nDistance to gene (kb)", 
                     breaks = c(1, 5, 8, 10, 12, 14, 17, 21), 
                     labels = c(gap_breaks[c(2, 5, 8)], c("5\'", "3\'"), gap_breaks[c(12, 15, 18)]) ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c(0, 0.2, 0.4, 0.6), 
                     limits = c(-0.1, 0.7)) +
  scale_color_manual(name = "Intergenic", values = colors_9[7:6]) +
  ylab("Gene transcription log2FC correlation") +
  theme_setting + 
  theme(strip.text = element_blank(),
        axis.ticks = element_line(), 
        panel.grid.major = element_line(), 
        panel.spacing = unit(1, "lines"),
        legend.position = "top")

dat_mat <- reshape2::melt(data = as.data.frame(mat[, c(1, 6:9)]), id = "Pos")
dat_mat$Sample <- rep(c(rep("2i", 21), rep("mTORi", 21)), 2)
dat_mat$Direction <- factor(c(rep("Enhancer", 21*2), rep("Other", 21*2)), 
                            levels = c("Enhancer", "Other"))
dat_mat$line_group <- factor(cumsum(is.na(dat_mat$value)))

g2 <- ggplot(dat_mat, aes(x = Pos, y = value, color = Direction, group = line_group)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey85") +
  # geom_vline(xintercept = 10, col = "grey85") +
  geom_rect(aes(xmin = 10, xmax = 12, ymin = -Inf, ymax = 0),
            fill = "#0000B2", linetype = 0) +
  geom_rect(aes(xmin = 9, xmax = 13, ymin = 0, ymax = Inf),
            fill = "grey95", linetype = 0) +
  geom_rect(aes(xmin = 10, xmax = 12, ymin = 0, ymax = Inf),
            fill = "grey75", linetype = 0) +
  geom_point() + 
  # geom_smooth(method = "loess", span = 0.45) +
  geom_line() +
  facet_grid(Sample~.) + 
  scale_x_continuous(name = "\nDistance to gene (kb)", 
                     breaks = c(1, 5, 8, 10, 12, 14, 17, 21), 
                     labels = c(gap_breaks[c(2, 5, 8)], c("5\'", "3\'"), gap_breaks[c(12, 15, 18)]) ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6),
                     labels = c(0, 0.2, 0.4, 0.6), 
                     limits = c(-0.1, 0.7)) +
  scale_color_manual(name = "", values = colors_9[c(4, 8)]) +
  ylab("") +
  theme_setting + 
  theme(strip.text = element_text(size=14),
        axis.ticks = element_line(),
        panel.grid.major = element_line(), 
        panel.spacing = unit(1, "lines"),
        legend.position = "top")

ggsave(plot = grid.arrange(g1, g2, ncol=2, widths = c(5, 5.5)),
       filename = "Fig2_TU_gene_log2FC_position_correlation_direction.png",
       path = "../fig2/figs", device = "png", width = 7, height = 5 )


# -----------------------------------------------------------------------------------------------------
# gene positioning co-expression test
# neighbor_gene_list <- list("downstream" = dse_genes, "divergent" = div_genes, "convergent" = con_genes)

TAD.gr <- importRanges('../data/TAD.Ren.total.HindIII.combined.domain_mm9.bed')
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(gene.gr$gene_id),
                       keytype = "ENTREZID",
                       columns = c("GENEID", "GENEBIOTYPE"))
gene.gr$gene_biotype <- res$GENEBIOTYPE[match(gene.gr$gene_id, res$ENTREZID)]
gene.gr$gene_id <- res$GENEID[match(gene.gr$gene_id, res$ENTREZID)]
gene.gr <- gene.gr[!is.na(gene.gr$gene_id) & !is.na(gene.gr$gene_biotype)]
gene.gr <- gene.gr[gene.gr$gene_biotype == 'protein_coding']
gene.gr <- gene.gr[width(gene.gr) < 2500000]

tad.mtch <- findOverlaps(TAD.gr, gene.gr, minoverlap = 200L)

tad.genes <- split(gene.gr$gene_id[subjectHits(tad.mtch)], queryHits(tad.mtch))

# read scRNA data
pool_cmp <- read.table('../../../sc_cell_cycle/data/mES_cmp_pool_count_log.txt', sep = '\t', header = T)

sc_exp_sl <- pool_cmp[, grep("_lif_", colnames(pool_cmp))]
sc_exp_sl <- sc_exp_sl[rownames(sc_exp_sl) %in% gene.gr$gene_id, ]
sc_exp_sl <- (t(exp(sc_exp_sl) - 1) / SizeFactorCal(exp(sc_exp_sl) - 1)) %>% log1p %>% t

sc_exp_2i <- pool_cmp[, grep("_2i_", colnames(pool_cmp))]
sc_exp_2i <- sc_exp_2i[rownames(sc_exp_2i) %in% gene.gr$gene_id, ]
sc_exp_2i <- (t(exp(sc_exp_2i) - 1) / SizeFactorCal(exp(sc_exp_2i) - 1)) %>% log1p %>% t
sc_exp_2i <- sc_exp_2i[, !is.na(colSums(sc_exp_2i))]

# covariance matrix
cov_mat_sl <- stats::cov(t(sc_exp_sl))
cov_mat_2i <- stats::cov(t(sc_exp_2i))

# extract covariance: divergent, convergent, downstream sense, TAD pairs, random background
get_cov <- function(gene.list, tad.genes, cov_mat) {
  res_list <- lapply(gene.list, 
                     function (x) 
                       cov_mat[x[apply(x, 1, function(y) all(y %in% colnames(cov_mat))), ]])
  names(res_list) <- names(gene.list)
  
  tad.genes <- lapply(tad.genes, function(x) x[x %in% colnames(cov_mat)])
  tad.genes <- tad.genes[lapply(tad.genes, length) %>% unlist > 1]
  
  tad_cov <- lapply(tad.genes, function(x)
  {
    cov_mat[t(combn(x, 2))]
  }
  ) %>% unlist
  
  bg_cov <- cov_mat[cbind(sample(colnames(cov_mat), 4000), sample(colnames(cov_mat), 4000)) ]
  c(res_list, list("TAD" = tad_cov, "Background" = bg_cov))
}

cov_res_sl <- get_cov(list("Divergent" = div_genes, 
                           "Convergent" = con_genes, 
                           "Down. sns." = dse_genes), 
                      tad.genes, cov_mat_sl)
cov_res_2i <- get_cov(list("Divergent" = div_genes, 
                           "Convergent" = con_genes, 
                           "Down. sns." = dse_genes), 
                      tad.genes, cov_mat_2i)
cov_res_sl <- cov_res_sl[c("Background", "TAD", "Down. sns.", "Convergent", "Divergent")]

# plot
plot_sc_cov_density(cov_res_sl, cov_res_2i, c("SL", "2i"))

ggsave(filename = "FigS2_gene_neighbor_coexpression.png",
       path = "../figS2/figs",
       device = "png", width = 4, height = 4.5 )

# p-values compare to background
wilcox.test(cov_res_sl$Background, cov_res_sl$TAD) # 0.204
wilcox.test(cov_res_sl$Background, cov_res_sl$`Down. sns.`) # 0.2464
wilcox.test(cov_res_sl$Background, cov_res_sl$Convergent) # 0.00000000000000022
wilcox.test(cov_res_sl$Background, cov_res_sl$Divergent) # 0.00000000000000022

# plot function
plot_sc_cov_density <- function(cov_res_list1, cov_res_list2, Samples) {
  dat1 <- data.frame(Covariance = unlist(cov_res_list1),
                    Positioning = rep(names(cov_res_list1), lengths(cov_res_list1)) )
  dat2 <- data.frame(Covariance = unlist(cov_res_list2),
                     Positioning = rep(names(cov_res_list2), lengths(cov_res_list2)) )
  dat <- rbind(dat1, dat2)
  dat$Sample <- factor(c(rep(Samples[1], nrow(dat1)), rep(Samples[2], nrow(dat2))),
                       levels = Samples)
  dat$Positioning <- factor(dat$Positioning, levels = names(cov_res_list1))
  
  ggplot(dat, aes(x=Covariance, fill=Sample)) +
    geom_vline(xintercept = 0, color = 'grey75', lwd = 2) +
    geom_histogram(aes(y=..density..), binwidth=.1, alpha=.7, position="identity") +
    geom_vline(data=plyr::ddply(dat, c("Positioning", "Sample"), summarise,
                                medians=median(Covariance)),
               aes(xintercept=medians, linetype = Sample),
               size=0.5) +
    facet_grid(Positioning~.) +
    scale_fill_manual(values = colors_20[c(13, 2)]) +
    xlim(c(-1, 2)) +
    ylab("Density") +
    theme_setting +
    theme(axis.text=element_text(size=10, face = "plain"),
          axis.title=element_text(size=14,face="bold"),
          strip.text.y = element_text(size = 12, face="bold")) 
}

# -------------------------------------------------------------------------------------------------
txRPK <- readRDS("../data/txRPK_SL_2i.RData")
tuRPK <- readRDS("../data/tuRPK_SL_2i.RData")

get_L_F_cor <- function(RPK) {
  L_mat <- RPK[, grep("LRNA", colnames(RPK))]
  L_mat <- sapply(unique(gsub("LRNA_(.*)_rep.", "\\1", colnames(L_mat))),
                  function(x) {
                    idx <- grep(x, colnames(L_mat))
                    .ifelse (length(idx) > 1, rowMeans(L_mat[, idx]), L_mat[, idx])
                  })
  F_mat <- RPK[, grep("FRNA", colnames(RPK))]
  F_mat <- sapply(unique(gsub("FRNA_(.*)_rep.", "\\1", colnames(F_mat))),
                  function(x) {
                    idx <- grep(x, colnames(F_mat))
                    .ifelse (length(idx) > 1, rowMeans(F_mat[, idx]), F_mat[, idx])
                  })
  keep.idx <- rowMeans(L_mat) > 0.05 # keep active genes
  message(paste("Number of TUs:", sum(keep.idx)))
  L_mat <- L_mat[keep.idx, ] %>% log() # add pseudo-count 0.01
  F_mat <- F_mat[keep.idx, ] %>% log()
  
  sapply(colnames(L_mat), function(x) single_variance_explained(L_mat[, x], F_mat[, x], T))
}

dat <- rbind(reshape::melt(get_L_F_cor(txRPK)), 
             reshape::melt(get_L_F_cor(tuRPK[grep("intergenic", rownames(tuRPK)), ])),
             reshape::melt(get_L_F_cor(tuRPK[grep("antisense", rownames(tuRPK)), ])),
             reshape::melt(get_L_F_cor(tuRPK[grep("uaRNA", rownames(tuRPK)), ])),
             reshape::melt(get_L_F_cor(tuRPK[grep("conRNA", rownames(tuRPK)), ])),
             reshape::melt(get_L_F_cor(tuRPK[grep("daRNA", rownames(tuRPK)), ])) )

dat$Type <- factor(rep(c("mRNA", "intergenic", "asRNA", "uaRNA", "conRNA", "daRNA"), each = 5),
                   levels = c("mRNA", "intergenic", "asRNA", "conRNA", "uaRNA", "daRNA"))
dat$Sample <- factor(rep(rownames(dat)[1:5], 6), levels = c("SL" , "2i_2d", "2i_7d", "mTORi_1d", "mTORi_2d") )

ggplot(dat[dat$Sample %in% c("SL" , "2i_2d", "mTORi_1d", "mTORi_2d"), ], 
       aes(x = Sample, y = value, fill = Type)) + 
  geom_bar(stat = "identity",position=position_dodge()) +
  ylim(c(0,1)) +
  xlab("") + ylab("log RPK correlation") + ggtitle("Nascent RNA ~ total RNA") +
  theme_setting + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "FigS2_LRNA_FRNA_TU_type_correlation.png",
       path = "../figS2/figs",
       device = "png", width = 5, height = 4 )


