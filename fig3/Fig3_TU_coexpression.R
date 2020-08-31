# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# -----------------------------------------------------------------------------------------------------
# load references
TU.DE.mm10.gr <- readRDS("../fig1/data/TU.DE.mm10.RData")
TU.DE.mm10.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$baseMean > 0 & 
                                 TU.DE.mm10.gr$location %in% c('protein_coding', 'intergenic')]

gene.gr <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
gene.gr <- `seqlevelsStyle<-`(gene.gr, "UCSC")
gene.gr <- gene.gr[gene.gr$gene_biotype == 'protein_coding']
gene.gr <-  gene.gr[width(gene.gr) < 2500000]

mtch <- findOverlaps(TU.DE.mm10.gr, gene.gr)
TU.DE.mm10.gr$gene_id <- NA
TU.DE.mm10.gr$gene_id[queryHits(mtch)] <- gene.gr$gene_id[subjectHits(mtch)]

# -----------------------------------------------------------------------------------------------------
# gene-gene relative position
TU.coding.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$location == "protein_coding"]
TU.coding.promoters <- promoters(TU.coding.gr, upstream = 2000, downstream = 0)
levels(strand(TU.coding.promoters)) <- c("-", "+", "*")
TU.coding.ds <- flank(TU.coding.gr, width = 2000, start = F)
levels(strand(TU.coding.ds)) <- c("-", "+", "*")

# get matches
ds_mtch <- findOverlaps(TU.coding.gr, promoters(TU.coding.gr, upstream = 5000, downstream = 0))
ua_mtch <- findOverlaps(TU.coding.gr, TU.coding.promoters)
con_mtch <- findOverlaps(TU.coding.gr, TU.coding.ds)

# downstream sense
q.hits <- queryHits(ds_mtch)
s.hits <- subjectHits(ds_mtch)

multi <- ds_mtch[q.hits %in% q.hits[duplicated(q.hits)] | s.hits %in% s.hits[duplicated(s.hits)] ]
multi <- multi[abs(queryHits(multi) - subjectHits(multi)) == 1, ]

uniq <- ds_mtch[!q.hits %in% q.hits[duplicated(q.hits)] & !s.hits %in% s.hits[duplicated(s.hits)] ]
dse_genes <- cbind("Upstream" = TU.coding.gr$gene_id[c(queryHits(multi), queryHits(uniq))],
                   "Downstream" = TU.coding.gr$gene_id[c(subjectHits(multi), subjectHits(uniq))])

# divergent
p.embeded <- subjectHits(ua_mtch) %in%
  (findOverlaps(TU.coding.promoters, TU.coding.gr, type = "within") %>% queryHits)
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

# -----------------------------------------------------------------------------------------------------
# gene-intergenic TU neighboring
TU.intergenic.gr <- TU.DE.mm10.gr[TU.DE.mm10.gr$location == "intergenic"]
TU.intergenic.gr$id <- seq_along(TU.intergenic.gr)
TU.coding.promoters <- promoters(TU.coding.gr, upstream = 2000, downstream = 0)
TU.coding.ds <- flank(TU.coding.gr, width = 2000, start = F)

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
    
    .ifelse = function(boo, e1, e2) if (boo) e1 else e2
    
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

gaps <- start(TU.intergenic.gr[gene.intergenic.pairs$id]) - 
  cbind(start(TU.coding.gr[gene.intergenic.pairs$gene_id]),
        end(TU.coding.gr[gene.intergenic.pairs$gene_id]))
idx <- ifelse(gene.intergenic.pairs$strand == "+", 
              ifelse(grepl("up", gene.intergenic.pairs$type), 1, 2), 
              ifelse(grepl("up", gene.intergenic.pairs$type), 2, 1))
gene.intergenic.pairs$gap <- 
  sapply(seq_len(nrow(gaps)), 
         function(x) gaps[x, idx[x]]
  ) * ifelse(gene.intergenic.pairs$strand == "+", 1, -1)

names(TU.coding.gr) <- TU.coding.gr$gene_id
tmp_table <- cbind(mcols(TU.coding.gr[gene.intergenic.pairs$gene_id])[, 2:6], 
                   mcols(TU.intergenic.gr[gene.intergenic.pairs$id])[, 2:6])
colnames(tmp_table) <- paste0(c(rep("Gene_", 5), rep("TU_", 5)), colnames(tmp_table))
gene.intergenic.pairs <- cbind(gene.intergenic.pairs, tmp_table)
rm(tmp_table)

# ---------------------------------------------------------------------
# find states of intergenic TU
ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub("(.*?)_(.*)", "\\2", ChromHMM$name)
TU.intergenic.gr$chromHMM <- NA
TU.intergenic.prmtr <- promoters(TU.intergenic.gr, 1, 0)
mtch <- findOverlaps(TU.intergenic.prmtr, ChromHMM)
TU.intergenic.gr$chromHMM[queryHits(mtch)] <- ChromHMM$name[subjectHits(mtch)]
gene.intergenic.pairs$TU_ChromHMM <- TU.intergenic.gr$chromHMM[gene.intergenic.pairs$id]

levels(strand(TU.intergenic.prmtr)) <- c("-", "+", "*")
bi_idx <- findOverlaps(TU.intergenic.gr, promoters(TU.intergenic.prmtr, upstream = 1000, downstream = 500)) %>% 
  countQueryHits() > 0
TU.intergenic.gr$Direction <- ifelse(bi_idx, "Bidirection", "Unidirection")
gene.intergenic.pairs$TU_direction <- TU.intergenic.gr$Direction[gene.intergenic.pairs$id]

# show co-expression
# gene.intergenic.pairs <- readRDS("data/gene.intergenic.pairs.RData")
gene.intergenic.pairs <- gene.intergenic.pairs[!apply(gene.intergenic.pairs[, c(7, 12)], 1, function(x) any(is.na(x))), ]

gap_breaks <- c(min(gene.intergenic.pairs$gap), 
                c(-400, -200, -100, -75, -50, -30, -20, -10, 0, 10, 20, 30, 50, 75, 100, 200, 400) * 1000,
                max(gene.intergenic.pairs$gap))

gene.intergenic.pairs$gap_class <- cut(gene.intergenic.pairs$gap, 
                                       breaks = gap_breaks, 
                                       labels = 1:18) %>% as.numeric()
gene.intergenic.pairs <- gene.intergenic.pairs[!is.na(gene.intergenic.pairs$gap_class), ]
# saveRDS(gene.intergenic.pairs, "data/gene.intergenic.pairs.RData")

# plot TU coverage
tmp.intergenic.pairs <- gene.intergenic.pairs[abs(gene.intergenic.pairs$gap) < 1e5, ]
tmp.gr <- TU.intergenic.gr[tmp.intergenic.pairs$id]
as.idx <- grepl("antisense", tmp.intergenic.pairs$type)

tmp.cov.gr <- IRanges(start = tmp.intergenic.pairs$gap + 1e5 - ifelse(as.idx, width(tmp.gr), 0), 
                      width = width(tmp.gr))

s.cov <- coverage(tmp.cov.gr[!as.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)
as.cov <- coverage(tmp.cov.gr[as.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)

d.idx <- grepl("Bidirection", tmp.intergenic.pairs$TU_direction)
s.idx <- tmp.intergenic.pairs$TU_ChromHMM %in% c("Enhancer", "ActivePromoter", "BivalentChromatin",
                                                 "StrongEnhancer", "WeakEnhancer",
                                                 "TranscriptionTransition", "TranscriptionElongation")

b.a.cov <- coverage(tmp.cov.gr[d.idx & s.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)
b.i.cov <- coverage(tmp.cov.gr[d.idx & !s.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)
u.a.cov <- coverage(tmp.cov.gr[!d.idx & s.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)
u.i.cov <- coverage(tmp.cov.gr[!d.idx & !s.idx], width = 2e5) %>% spline(x = seq_len(2e5), n = 200) %>% "$"(y)

insert_NA <- function(x, p = 100, n_0 = 30) c(x[seq_len(p)], rep(NA, n_0), x[seq_len(p) + p])
dat <- data.frame(Pos = rep(1:230, 6),
                  Cov = c(b.a.cov %>% insert_NA(), 
                          b.i.cov %>% insert_NA(),
                          u.a.cov %>% insert_NA(),
                          u.i.cov %>% insert_NA(),
                          s.cov %>% insert_NA(),
                          as.cov %>% insert_NA()),
                  Type = c(rep("Bidirectional active", 230), rep("Bidirectional inactive", 230),
                           rep("Unidirectional active", 230), rep("Unidirectional inactive", 230),
                           rep("Sense", 230), rep("Antisense", 230)))
plot_tu_cov <- function(dat, type, l_color) {
  ggplot(dat[dat$Type == type, ], aes(x = Pos, y = Cov)) +
    geom_rect(aes(xmin = 100, xmax = 130, ymin = -Inf, ymax = 0),
              fill = "#0000B2", linetype = 0) +
    geom_rect(aes(xmin = 100, xmax = 130, ymin = 0, ymax = Inf),
              fill = "grey", linetype = 0) +
    geom_line(lwd = 2, color = l_color) +
    scale_x_continuous(name = "\nDistance from gene (kb)\n",
                       breaks = c(0, 50, 100, 130, 180, 230), 
                       labels = c("-100", "-50", "Start", "End", "50", "100")) +
    ylab("Occurence") +
    ggtitle(type) +
    theme_setting
}

g1 <- plot_tu_cov(dat, type = "Bidirectional active", l_color = colors_20[1])
g2 <- plot_tu_cov(dat, type = "Bidirectional inactive", l_color = colors_20[6])
g3 <- plot_tu_cov(dat, type = "Unidirectional active", l_color = colors_20[15])
g4 <- plot_tu_cov(dat, type = "Unidirectional inactive", l_color = colors_20[5])
g5 <- plot_tu_cov(dat, type = "Sense", l_color = colors_20[8])
g6 <- plot_tu_cov(dat, type = "Antisense", l_color = colors_20[10])

ggsave(plot = grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 2),
       filename = "FigS3_TU_gene_neighbor_occurence_direction_state.png", path = "../figS3/figs",
       device = "png", width = 8.5, height = 9 )

# plot by TU chromatin state
mat <- data.frame()
for( i in 1:18 ) {
  tmp.pairs <- gene.intergenic.pairs[as.numeric(gene.intergenic.pairs$gap_class) == i, ]
  s.idx <- tmp.pairs$TU_ChromHMM %in% c("Enhancer", "ActivePromoter", 
                                                    "StrongEnhancer", "WeakEnhancer")
  d.idx <- tmp.pairs$TU_direction == "Bidirection"
  tmp.pairs <- tmp.pairs[, grep("log2", colnames(gene.intergenic.pairs))] %>% as.matrix()
  mat <- rbind(mat, c(i, 
                      cor(tmp.pairs[s.idx & d.idx, ])[rbind(c(1,3), c(2,4))],
                      cor(tmp.pairs[s.idx & !d.idx, ])[rbind(c(1,3), c(2,4))],
                      cor(tmp.pairs[!s.idx & d.idx, ])[rbind(c(1,3), c(2,4))],
                      cor(tmp.pairs[!s.idx & !d.idx, ])[rbind(c(1,3), c(2,4))]
                      ))
}
# colnames(mat) <- c("Pos", "Cor_Inactive_2i", "Cor_Active_2i", "Cor_Inactive_mTORi", "Cor_Active_mTORi")
mat <- reshape2::melt(data = mat, id = "X1")
mat$Sample <- rep(c(rep("2i", 18), rep("mTORi", 18)), 4)
mat$Direction <- rep(c(rep("Unidirection", 18*2), rep("Bidirection", 18*2)), 2)
mat$State <- c(rep("Active", 18*4), rep("Inactive", 18*4))

ggplot(mat, aes(x = X1, y = value, color = Direction)) +
  geom_point() + 
  geom_smooth(method = "loess", span = 0.45) +
  geom_hline(yintercept = 0, col = "grey") +
  facet_grid(Sample~State) + 
  scale_x_continuous(name = "Distance from gene (kb)", breaks = c(1, 5, 8, 10, 12, 15, 18), 
                     labels = gap_breaks[c(2, 5, 8, 10, 12, 15, 18)] / 1000) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(-0.4, 1)) +
  scale_color_manual(name = "Intergenic TU", values = colors_9[7:6]) +
  ylab("TU~gene Log2FC correlation") +
  ggtitle("Neighborhood coexpression") +
  theme_setting + 
  theme(strip.text = element_text(size=14),
        axis.ticks = element_line(),
        panel.spacing = unit(1, "lines"))
ggsave(filename = "Fig3_TU_gene_log2FC_position_correlation_direction_state.png",
       path = "../fig3/figs", device = "png", width = 8.5, height = 6 )

# ---------------------------------------------------------------------
sample_ncRNA_counts_Rates <- readRDS("../figS2/data/sample_nc_counts_Rates_combined.RData")
tu_names <- paste("intergenic", seqnames(TU.intergenic.gr), 
                  start(TU.intergenic.gr), end(TU.intergenic.gr), 
                  strand(TU.intergenic.gr), sep = "_")
intergenic_sl_count_rates <- sample_ncRNA_counts_Rates %>%
  dplyr::filter(Sample == "SL") %>%
  dplyr::filter(gene_id %in% tu_names) 
intergenic_2i_count_rates <- sample_ncRNA_counts_Rates %>%
  dplyr::filter(Sample == "2i_2d") %>%
  dplyr::filter(gene_id %in% tu_names) 
intergenic_mTORi_count_rates <- sample_ncRNA_counts_Rates %>%
  dplyr::filter(Sample == "mTORi_1d") %>%
  dplyr::filter(gene_id %in% tu_names) 

TU.intergenic.gr$half_life_2i <- 
  intergenic_2i_count_rates$Half_life[match(tu_names, intergenic_sl_count_rates$gene_id)] %>% 
  log2()
TU.intergenic.gr$half_life_mTORi <- 
  intergenic_mTORi_count_rates$Half_life[match(tu_names, intergenic_sl_count_rates$gene_id)] %>% 
  log2()

TU.intergenic.gr$half_life_delta_2i <- TU.intergenic.gr$half_life %>% 
  "/"(intergenic_2i_count_rates$Half_life[match(tu_names, intergenic_2i_count_rates$gene_id)]) %>%
  log2()
TU.intergenic.gr$half_life_delta_mTORi <- TU.intergenic.gr$half_life %>% 
  "/"(intergenic_mTORi_count_rates$Half_life[match(tu_names, intergenic_mTORi_count_rates$gene_id)]) %>%
  log2()

gene.intergenic.pairs <-
  cbind(gene.intergenic.pairs, mcols(TU.intergenic.gr)[gene.intergenic.pairs$id,
                                                       c("half_life_2i",
                                                         "half_life_mTORi",
                                                         "half_life_delta_2i", 
                                                         "half_life_delta_mTORi")])
gene.intergenic.pairs <- gene.intergenic.pairs[complete.cases(gene.intergenic.pairs), ]

mat <- data.frame(matrix(ncol = 9, nrow = 0))
for( i in 1:18 ) {
  tmp.pairs <- gene.intergenic.pairs[as.numeric(gene.intergenic.pairs$gap_class) == i, ]
  # s.idx <- tmp.pairs$TU_ChromHMM %in% c("Enhancer", "ActivePromoter", 
  #                                       "StrongEnhancer", "WeakEnhancer")
  s.idx <- tmp.pairs$TU_direction == "Bidirection"
  tmp.pairs <- tmp.pairs[, c(7, 9, 19:22)] %>% as.matrix()
  res_2i_s_T <- cor(tmp.pairs[s.idx, -c(4, 6)])
  res_2i_s_F <- cor(tmp.pairs[!s.idx, -c(4, 6)])
  res_mTORi_s_T <- cor(tmp.pairs[s.idx, -c(3, 5)])
  res_mTORi_s_F <- cor(tmp.pairs[!s.idx, -c(3, 5)])
  
  mat <- rbind(mat, c(i, 
                      res_2i_s_T[1,3], res_2i_s_F[1,3],
                      res_2i_s_T[2,4], res_2i_s_F[2,4],
                      res_mTORi_s_T[1,3], res_mTORi_s_F[1,3], 
                      res_mTORi_s_T[2,4], res_mTORi_s_F[2,4]
  ))
}
colnames(mat) <- c("Pos", 
                   "2i|Half-life|Bidirection",
                   "2i|Half-life|Unidirection",
                   "2i|Half-life change|Bidirection",
                   "2i|Half-life change|Unidirection",
                   "mTORi|Half-life|Bidirection",
                   "mTORi|Half-life|Unidirection",
                   "mTORi|Half-life change|Bidirection",
                   "mTORi|Half-life change|Unidirection")
mat <- reshape2::melt(data = mat, id = "Pos")
mat <- cbind(mat, stringr::str_split_fixed(mat$variable, "\\|", 3))
colnames(mat)[4:6] <- c("Sample", "Type", "State")

ggplot(mat, aes(x = Pos, y = value, color = State)) +
  geom_point() + 
  geom_smooth(method = "loess", span = 0.45) +
  geom_hline(yintercept = 0, col = "grey") +
  facet_grid(Sample~Type) + 
  scale_x_continuous(name = "Distance from gene (kb)", breaks = c(1, 5, 8, 10, 12, 15, 18), 
                     labels = gap_breaks[c(2, 5, 8, 10, 12, 15, 18)] / 1000) +
  scale_y_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), labels = c(-0.6, -0.3, 0, 0.3, 0.6), limits = c(-0.6, 0.6)) +
  scale_color_manual(name = "Intergenic TU", values = colors_9[c(8, 2)]) +
  ylab("TU half-life ~gene Log2FC correlation") +
  ggtitle("Neighborhood coexpression") +
  theme_setting + 
  theme(strip.text = element_text(size=14),
        axis.ticks = element_line(),
        panel.spacing = unit(1, "lines"))

ggsave(filename = "Fig3_TU_half_life_gene_log2FC_position_correlation_direction.png",
       path = "../fig3/figs", device = "png", width = 8.5, height = 6 )

# ---------------------------------------------------------------------
if (FALSE) {
  # test if intergenic TU links gene coexpression 
  intermediate_TU_id <- gene.intergenic.pairs$id[duplicated(gene.intergenic.pairs$id) & abs(gene.intergenic.pairs$gap) < 100000]
  gene.intergenic.gene.pairs <- gene.intergenic.pairs[gene.intergenic.pairs$id %in% intermediate_TU_id, ]
  gene.intergenic.gene.pairs <- gene.intergenic.gene.pairs[order(gene.intergenic.gene.pairs$id), ]
  
  all_sided <- function(x, r) all(x < r[1]) | all(x > r[2])
  link_dat <- foreach(i = unique(gene.intergenic.gene.pairs$id), .combine = rbind) %dopar% {
    tmp_table <- gene.intergenic.gene.pairs[gene.intergenic.gene.pairs$id == i, ]
    
    c(i, sum(abs(tmp_table$Gene_log2FoldChange_2i - tmp_table$TU_log2FoldChange_2i)),
      sum(abs(tmp_table$Gene_log2FoldChange_mTORi - tmp_table$TU_log2FoldChange_mTORi)),
      sd(c(tmp_table$Gene_log2FoldChange_2i, tmp_table$TU_log2FoldChange_2i[1])),
      sd(c(tmp_table$Gene_log2FoldChange_mTORi, tmp_table$TU_log2FoldChange_mTORi[1])),
      c(tmp_table$Gene_log2FoldChange_2i, tmp_table$TU_log2FoldChange_2i[1]) %>% all_sided(c(-0.7, 0.7)),
      c(tmp_table$Gene_log2FoldChange_mTORi, tmp_table$TU_log2FoldChange_mTORi[1]) %>% all_sided(c(-0.7, 0.7)) )
  }
  
  link_dat <- data.frame(link_dat)
  colnames(link_dat) <- c("TU_id", "Exp_diff_2i", "Exp_diff_mTORi", "Exp_sd_2i", "Exp_sd_mTORi", "Same_diff_2i", "Same_diff_mTORi")
  link_dat$Gene_diff_2i <- aggregate(gene.intergenic.gene.pairs$Gene_log2FoldChange_2i,
                                     list(gene.intergenic.gene.pairs$id), 
                                     function(x) abs(diff(x)))$x
  link_dat$Gene_diff_mTORi <- aggregate(gene.intergenic.gene.pairs$Gene_log2FoldChange_mTORi,
                                        list(gene.intergenic.gene.pairs$id), 
                                        function(x) abs(diff(x)))$x
  link_dat$Gene_mean_2i <- aggregate(gene.intergenic.gene.pairs$Gene_log2FoldChange_2i,
                                     list(gene.intergenic.gene.pairs$id), 
                                     mean)$x
  link_dat$Gene_mean_mTORi <- aggregate(gene.intergenic.gene.pairs$Gene_log2FoldChange_mTORi,
                                        list(gene.intergenic.gene.pairs$id), 
                                        mean)$x
  
  link_dat <- cbind(link_dat, t(matrix(gene.intergenic.gene.pairs$gene_id, nrow = 2)))
  link_dat$Direction <- TU.intergenic.gr$Direction[link_dat$TU_id]
  
  table(link_dat$Direction, link_dat$Same_diff_2i) 
  
  intergenic_linked_genes <- cbind(as.character(link_dat$`1`[link_dat$Same_diff_2i == 1]),
                                   as.character(link_dat$`2`[link_dat$Same_diff_2i == 1])) %>% 
    apply(1, base::sort) %>% t() %>% as.data.frame()
  intergenic_linked_genes <- intergenic_linked_genes[!duplicated(intergenic_linked_genes), ]
  
  intergenic_not_linked_genes <- cbind(as.character(link_dat$`1`[link_dat$Same_diff_2i == 0]),
                                       as.character(link_dat$`2`[link_dat$Same_diff_2i == 0])) %>% 
    apply(1, base::sort) %>% t() %>% as.data.frame()
  intergenic_not_linked_genes <- intergenic_not_linked_genes[!duplicated(intergenic_not_linked_genes), ]
}
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
                           "Down. cis" = dse_genes), 
                      tad.genes, cov_mat_sl)
cov_res_2i <- get_cov(list("Divergent" = div_genes, 
                           "Convergent" = con_genes, 
                           "Down. cis" = dse_genes), 
                      tad.genes, cov_mat_2i)
cov_res_sl <- cov_res_sl[c("Background", "TAD", "Down. cis", "Convergent", "Divergent")]
# plot
plot_sc_cov_density(cov_res_sl, cov_res_2i, c("SL", "2i"))
ggsave(filename = "FigS3_gene_neighbor_coexpression.png", path = "../figS3/figs",
       device = "png", width = 4, height = 6 )

# p-values compare to background
wilcox.test(cov_res_sl$Background, cov_res_sl$TAD) # 0.204
wilcox.test(cov_res_sl$Background, cov_res_sl$`Down. cis`) # 0.2464
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
    geom_histogram(aes(y=..density..), binwidth=.1, alpha=.7, position="identity") +
    geom_vline(data=plyr::ddply(dat, c("Positioning", "Sample"), summarise,
                                medians=median(Covariance)),
               aes(xintercept=medians),
               linetype = 2, size=0.5) +
    facet_grid(Positioning~.) +
    # scale_fill_manual(values = c_col) +
    xlim(c(-1, 2)) +
    ylab("Density") +
    theme_setting +
    theme(axis.text=element_text(size=10, face = "plain"),
          axis.title=element_text(size=14,face="bold"),
          strip.text.y = element_text(size = 12, face="bold")) 
}

# ----------------------------------------------------------------------------------------
rho = readRDS("../fig2/data/est.rho.RData")
high_genes = names(tail(sort(rho), 400))
low_genes = names(head(sort(rho), 400))

high_genes <- high_genes[high_genes %in% rownames(sc_exp_sl)]
low_genes <- low_genes[low_genes %in% rownames(sc_exp_sl)]

cv_h = sapply(high_genes, function(x) sd(exp(sc_exp_sl[x, ])) / mean(exp(sc_exp_sl[x, ])))
cv_l = sapply(low_genes, function(x) sd(exp(sc_exp_sl[x, ])) / mean(exp(sc_exp_sl[x, ])))

dat = data.frame(Saturation = c(rep("High", length(cv_h)), rep("Low", length(cv_l))),
                 Coefficient.Variation = c(cv_h, cv_l))

ggplot(dat, aes(x = Coefficient.Variation, fill = Saturation)) +
  geom_histogram( binwidth = 0.5, alpha=.7, position="identity") +
  geom_vline(data=plyr::ddply(dat, "Saturation", summarise,
                              medians=median(Coefficient.Variation)),
             aes(xintercept = medians),
             linetype = 2, size= 1, col = colors_20[c(1,8)]) +
  geom_text(x = 14, y = 50, label = paste("p <", round(wilcox.test(cv_h, cv_l)$p.val, 2)) ) +
  scale_fill_manual(values = colors_20[c(1,8)]) +
  theme_setting
ggsave(filename = "FigS2_scRNA_coefficient_variation_rho_high_low.png", path = "../figS2/figs",
       device = "png", width = 5, height = 3 )
