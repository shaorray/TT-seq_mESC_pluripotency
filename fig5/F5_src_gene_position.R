# get gene types: convergent, divergent, downstream

require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

ref.genes <- genes(txdb)[findOverlaps(genes(txdb), gene.gr, ignore.strand = F) %>%
                           countQueryHits() > 0] %>% sort


res <- res[!grepl("Mir|Snord35", res$SYMBOL) & !duplicated(res$SYMBOL), ]
ref.genes <- ref.genes[ref.genes$gene_id %in% res$ENTREZID]

ref.genes$gene_id = res$ENSEMBL[match(ref.genes$gene_id, res$ENTREZID)]
ref.genes <- ref.genes[findOverlaps(ref.genes, ref.genes, type = "within") %>%
                         countQueryHits() == 1]
ref.genes <- ref.genes[width(ref.genes) < 1000000]

con_genes <- NULL # convergent
div_genes <- NULL # divergent
dse_genes <- NULL # downstream sense

# downstream sense
.promoters <- promoters(ref.genes, upstream = 5000, downstream = 0)
mtch <- findOverlaps(ref.genes, .promoters)

q.hits <- queryHits(mtch)
s.hits <- subjectHits(mtch)

multi <- mtch[q.hits %in% q.hits[duplicated(q.hits)] | s.hits %in% s.hits[duplicated(s.hits)] ]
multi <- multi[abs(queryHits(multi) - subjectHits(multi)) == 1, ]

uniq <- mtch[!q.hits %in% q.hits[duplicated(q.hits)] & !s.hits %in% s.hits[duplicated(s.hits)] ]
dse_genes <- cbind("Upstream" = ref.genes$gene_id[c(queryHits(multi), queryHits(uniq))],
                   "Downstream" = ref.genes$gene_id[c(subjectHits(multi), subjectHits(uniq))])

# divergent
.promoters <- promoters(ref.genes, upstream = 2000, downstream = 0)
levels(strand(.promoters)) <- c("-", "+", "*")
mtch <- findOverlaps(ref.genes, .promoters)

p.embeded <- subjectHits(mtch) %in%
  (findOverlaps(.promoters, ref.genes, type = "within") %>% queryHits)
mtch <- mtch[!p.embeded]

q.hits <- queryHits(mtch)
s.hits <- subjectHits(mtch)
multi <- mtch[q.hits %in% q.hits[duplicated(q.hits)] | s.hits %in% s.hits[duplicated(s.hits)] ]
ms.hits <- subjectHits(multi)
p.multi <- ms.hits[duplicated(ms.hits)]
multi.mtch <- multi[ms.hits %in% p.multi]
multi.strd <- ref.genes[p.multi] %>% strand %>% as.character

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

uniq <- mtch[!q.hits %in% q.hits[duplicated(q.hits)] & !s.hits %in% s.hits[duplicated(s.hits)] ]
div_genes <- cbind("Promoter" = ref.genes$gene_id[c(p.multi, queryHits(uniq))],
                   "Gene" = ref.genes$gene_id[c(p.closest, subjectHits(uniq))])

# convergent
.termination <- flank(ref.genes, width = 2000, start = F)

levels(strand(.termination)) <- c("-", "+", "*")
mtch <- findOverlaps(ref.genes, .termination)

q.hits <- queryHits(mtch)
s.hits <- subjectHits(mtch)
multi <- mtch[q.hits %in% q.hits[duplicated(q.hits)] ]

f.multi <- unique(queryHits(multi))
multi.strd <- ref.genes[f.multi] %>% strand %>% as.character

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

uniq <- mtch[ !q.hits %in% q.hits[duplicated(q.hits)] ]
con_genes <- cbind("Forward" = ref.genes$gene_id[c(f.multi, queryHits(uniq))],
                   "Reverse" = ref.genes$gene_id[c(f.closest, subjectHits(uniq))])