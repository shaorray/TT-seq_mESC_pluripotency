# get last exons of each gene

if (!file.exists("data/gene.last.exon.mm9.gtf"))
{
  # get the last exons
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
  seqlevels(txdb) <- paste0("chr", c(1:19, "X", "Y")) 
  
  # processing refers to G4hunter
  transcripts <- transcripts(txdb)	#55105
  names(transcripts) <- transcripts$tx_name
  transBgen <- transcriptsBy(txdb, by="gene")	#23434 more than gene number!
  transBgenClean <- transBgen[names(transBgen) %in% names(genes(txdb))]	#23274
  transclean <- transcripts[transcripts$tx_name %in% unlist(transBgenClean)$tx_name]
  
  txnames <- transclean$tx_name
  
  exonallnames <- exonsBy(txdb,by='tx', use.name=T)	#55105
  exonallclean <- exonallnames[names(exonallnames) %in% txnames]	#47847
  
  exonlast <- unlist(exonallclean)[cumsum(lengths(exonallclean))] #47847
  
  # remove redundant and silent exons
  gene.gr <- readRDS("../fig4/data/gene.gr.RData")
  gene.exonlast <- exonlast[findOverlaps(exonlast, gene.gr, ignore.strand = F) %>% countQueryHits() > 0]
  
  # convert UCSC tx id to emsembl gene id
  res.txdb <- biomaRt::select(TxDb.Mmusculus.UCSC.mm9.knownGene,
                              keys = names(uni.exonlast),
                              keytype = "TXNAME",
                              columns = "GENEID")
  library(org.Mm.eg.db)
  ks <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
  res <- biomaRt::select(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL",
                         columns = c("ENTREZID", "SYMBOL"))
  gene.exonlast$gene_id = res$ENSEMBL[match(res.txdb$GENEID[match(names(gene.exonlast), res.txdb$TXNAME)],
                                           res$ENTREZID)]
  # keep genes existing in annotation
  gene.exonlast <- gene.exonlast[gene.exonlast$gene_id %in% gsub("\\..*", "\\2", gene.gr$gene_id)]
  
  # save data
  export.gff3(gene.exonlast, "data/gene.last.exon.mm9.gff3")
  require(GenomicFeatures)
  gene.promoter.gr <- promoters(genes(txdb), upstream = 1, downstream = 0)
  terWindow <- foreach ( i = seq_along(gene.gr), .combine = c) %dopar%
  {
    tmp.terWindow <- flank(gene.gr[i], width = 15000, start = F)
    tmp.promoter.gr <- gene.promoter.gr[queryHits(findOverlaps(gene.promoter.gr, tmp.terWindow))]
    
    if (length(tmp.promoter.gr) > 0)
    {
      if (as.character(strand(tmp.terWindow)) == "+")
      {
        end(tmp.terWindow) <- start(tmp.promoter.gr)[which.min(start(tmp.promoter.gr) - start(tmp.terWindow))]
      } else {
        start(tmp.terWindow) <- start(tmp.promoter.gr)[which.min(end(tmp.terWindow) - start(tmp.promoter.gr))]
      }
    }
    tmp.terWindow
  }
  terWindow$gene_id <- gene.gr$gene_id
  export.gff3(terWindow, "data/terWindow.mm9.gff3")
}

