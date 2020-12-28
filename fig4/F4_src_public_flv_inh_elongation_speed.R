# Gro-seq samples from Jonkers et. al 2014
# TU annotation raw intervals (without joining by gene exon overlaps) 
# flv inhibition for 2, 5, 12.5, 25, 50 minutes

# use annotated TUs temporary objects from TU filter as the current tx 5' end positions -----------------------------------------
anno.files = list.files(path = "../data/TU_anno/Gro_seq_FP", pattern = "gff3", full.names = T)
all.TU.gr.list = sapply(anno.files, importRanges)
names(all.TU.gr.list) = lapply(anno.files, function(x) gsub('.*TU_filter\\+(.*)(minFP|).gff3', '\\1', x[1])) %>% unlist()
all.TU.gr.list = all.TU.gr.list[c(3, 5, 1, 2, 4)]

# get reference genes
gene_id_mm10 <- sapply(all.TU.gr.list, function(x) x$gene_id) %>% 
  unlist() %>% unique() %>% "["( complete.cases(.))

gene.mm10 <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
gene.mm10 <- gene.mm10[
  # gene.mm10$gene_id %in% gene_id_mm10 &
                         gene.mm10$gene_biotype %in% c('protein_coding', 'lncRNA') &
                         width(gene.mm10) > 3500]
gene.mm10 <- `seqlevelsStyle<-`(gene.mm10, "UCSC")

registerDoParallel(cores = 30)

dist_table <- foreach( i = seq_along(gene.mm10), .combine = rbind) %dopar% {
  
  gene.tmp <- gene.mm10[i]
  tss.tmp <- promoters(gene.tmp, upstream = 1000, downstream = 1000)
  gene.strand <- as.character(strand(gene.tmp))
  gene.tu.tmp <- lapply(all.TU.gr.list, 
                     function(x) x[findOverlaps(x, gene.tmp) %>% countQueryHits() > 0]) %>% 
    Reduce("c", .) %>% reduce() %>% sort()
  
  append_granges <- function(x, gene) {
    # return appended TTS boundary of gene intervals
    stopifnot(length(unique(seqnames(x))) == 1)
    x <- sort(x)
    end(x[1]) <- end(x[length(x)])
    gene.out <- gene
    if (as.character(strand(gene) == "+")) {
      if (end(gene.out) < end(x[length(x)]))
        end(gene.out) <- end(x[length(x)])
    } else {
      if (start(gene.out) > start(x[1]))
        start(gene.out) <- start(x[1])
    }
    gene.out
  }
  
  if (sum(width(gene.tu.tmp)) / width(gene.tmp) < 0.3 | # skip low coverage genes
       sum(findOverlaps(gene.tu.tmp, gene.mm10) %>% countSubjectHits() > 1)) {
    return(rep(0, length(all.TU.gr.list)))
  } else {
    gene.tu.tmp <- append_granges(gene.tu.tmp, gene.tmp)
  }
  
  gap <- NULL
  for ( j in seq_along(all.TU.gr.list) )
  {
    TU.tmp <- all.TU.gr.list[[j]]
    TU.tmp <- TU.tmp[findOverlaps(TU.tmp, gene.tu.tmp) %>% countQueryHits > 0]
    
    # TU.tmp <- TU.tmp[width(TU.tmp) > 2000]
    
    # exclude TSS downstream 2kb region for overlapping TUs, since GRO-seq TSS peaks are consistent
    # but except the TU that is the longest of all called intervals for a given gene
    if (length(TU.tmp) > 1) # remove the TSS peak if exists
    {
      TU.tmp <- TU.tmp[findOverlaps(TU.tmp, tss.tmp) %>% countQueryHits() == 0 |
                         (width(TU.tmp) / width(gene.tu.tmp) > 0.05)]
    }
    
    if (length(TU.tmp) > 1) # remove small truncated TUs
    {
      TU.tmp <- TU.tmp[findOverlaps(TU.tmp, TU.tmp + 1000) %>% countQueryHits() > 1 |
                         (width(TU.tmp) / width(gene.tu.tmp) > 0.05)]
    }
    
    if ( gene.strand == '+' ) 
    {
      tmp_gap <- start(TU.tmp) - start(gene.tu.tmp) 
    } else {
      tmp_gap <- end(gene.tu.tmp) - end(TU.tmp) 
    }
    
    gap <- c(gap, min(tmp_gap))
  }
  gap[gap < 0] <- 0
  gap
}

rownames(dist_table) <- gene.mm10$gene_id
dist_table <- dist_table - dist_table[, 1]
dist_table2 <- dist_table[!apply(dist_table, 1, 
                                function(x) all(is.infinite(x)) | is.na(x[1]) | sum(x) <= 0), ]
dist_table2[dist_table2 < 0] <- 0

# find pause-release genes
# sum_dist <- apply(dist_table, 1, function(x) x[!is.infinite(x) & x > (-10000)] %>% diff %>% sum)


# compute elongation speed with 5, 12.5, 25, 50 minutes -------------------------------------------------------------------
# intercept adjusts lag time of each individual gene in responses to flv inhibition

durations <- matrix(c(2, 5, 12.5, 25, 50), ncol = 1)
elongation_speed_table <- data.frame(speed = rep(0, nrow(dist_table2)), # minute / kb
                               lag_time = rep(0, nrow(dist_table2)), # minute
                               time_points = rep(0, nrow(dist_table2)), # number of time points used for speed estimation
                               row.names = rownames(dist_table2))

for ( i in seq_len(nrow(dist_table2)) )
{
  idx <- which(!is.infinite(dist_table2[i, ]) & dist_table2[i, ] > 0)
  if (dist_table2[i, 4] > dist_table2[i, 5]) idx <- idx[!idx == 5]
  if (length(unique(dist_table2[i, idx])) == 1) idx <- min(idx)
  
  elongation_speed_table$time_points[i] <- length(idx)
  
  if (length(idx) > 1) {
    x <- cbind(1, dist_table2[i, idx] / 1000) # change unit to kb
    betas <- solve(t(x) %*% x) %*% t(x) %*% durations[idx, ]
    if (betas[1, ] < 0) {
      elongation_speed_table$speed[i] <- mean(dist_table2[i, idx] / 1000 / durations[idx, ])
    } else {
      elongation_speed_table$lag_time[i] <- betas[1, ]
      elongation_speed_table$speed[i] <- 1 / betas[2, ]
    }
    
  } else if (length(idx) == 1) {
    elongation_speed_table$speed[i] <- dist_table2[i, idx] / 1000 / durations[idx, ]
  }
}

elongation_speed_table <- elongation_speed_table[elongation_speed_table$speed > 0, ] 

write.table(elongation_speed_table, '../fig4/data/elongation_speed_table.txt',
            sep = '\t', quote = F, row.names = T, col.names = T)
