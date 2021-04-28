library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggpubr)

enrichGeneSets <- function(gene_id, method = "GO", title = NULL, 
                           ontology = "CC", top_n_go = 10, colorset = viridis::viridis(10))
{ # Agrs:
  # gene_id: ENSEMBL gene id of interest
  # ontology: BP = Biological Process; CC = Cellular Component; MF = Molecular Function
  
  all_gene_ids <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)$gene_id
  
  
  
  # convert to ENTREZ_GENE_ID
  ks <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
  res <- biomaRt::select(org.Mm.eg.db, keys = ks, keytype = "ENSEMBL", 
                         columns = c("ENTREZID", "SYMBOL"))
  
  gene_id <- res$ENTREZID[match(gene_id, res$ENSEMBL)]
  
  if (method == 'GO')
  { # 
    ego <- enrichGO(gene          = gene_id,
                    universe      = all_gene_ids,
                    OrgDb         = org.Mm.eg.db,
                    ont           = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  }
  
  if (method == 'DAVID')
  {
    # if(!require(RDAVIDWebService)) BiocManager::install("RDAVIDWebService")
    # if(!require(rJava)) install.packages("rJava")
    # library(rJava)
    # library(RDAVIDWebService)
    # 
    # ego <- enrichDAVID(gene          = gene_id,
    #                    idType        = "ENTREZ_GENE_ID",
    #                    universe      = all_gene_ids,
    #                    minGSSize     = 10,
    #                    maxGSSize     = 500,
    #                    annotation    = "GOTERM_BP_FAT",
    #                    pvalueCutoff  = 0.05,
    #                    pAdjustMethod = "BH", 
    #                    qvalueCutoff  = 0.2)
  }
  
  
  if (method == 'enricher')
  {
    # ego <- enricher(gene             = gene_id, 
    #                 pvalueCutoff     = 0.05, 
    #                 pAdjustMethod    = "BH", 
    #                 universe         = all_gene_ids,
    #                 minGSSize        = 10, 
    #                 maxGSSize        = 500, 
    #                 qvalueCutoff     = 0.2, 
    #                 TERM2GENE,
    #                 TERM2NAME = NA)
  }
  
  if (method == 'KEGG')
  {
    ego <- enrichKEGG(gene              = gene_id, 
                      organism          = "mmu", 
                      keyType           = "kegg",
                      pvalueCutoff      = 0.05, 
                      pAdjustMethod     = "BH", 
                      universe          = all_gene_ids,
                      minGSSize         = 10, 
                      maxGSSize         = 500, 
                      qvalueCutoff      = 0.2,
                      use_internal_data = FALSE)
  }
  
  if (method == 'MKEGG')
  {
    ego <- enrichMKEGG(gene              = gene_id, 
                       organism          = "mmu", 
                       keyType           = "kegg",
                       pvalueCutoff      = 0.05, 
                       pAdjustMethod     = "BH", 
                       universe          = all_gene_ids,
                       minGSSize         = 10, 
                       maxGSSize         = 500, 
                       qvalueCutoff      = 0.2)
  }
  
  if (is.null(title)) title = method
  
  idx <- order(ego@result$Count, decreasing = T)[seq_len(top_n_go)]
  ego@result$pvalue
  dat <- with(ego@result[idx, ], 
              data.frame(Description = Description, 
                         GO = ID,
                         p.adjust = p.adjust, 
                         GeneRatio = Count / length(ego@gene),
                         Count = Count
              )
  )
  dat$p.adjust <- as.numeric(dat$p.adjust)
  dat$GeneRatio <- as.numeric(dat$GeneRatio)
  
  
  ggdotchart(dat, x = "Description", y = "GeneRatio",
             color = "p.adjust",                                 # Color by groups
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             rotate = T,
             # group = "cyl",                                # Order by groups
             dot.size = 7,                                 # Large dot size
             label = round(dat$Count),                        # Add counts as dot labels
             font.label = list(color = "grey", size = 10, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr()                        # ggplot2 theme
  ) + 
    scale_color_gradientn(colours = colorset, limits=c(0, 1)) + 
    theme_setting +
    theme(legend.position = "right", 
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = unit(c(1,1,1,1), "cm")) +
    xlab("") + 
    ggtitle(title)
     
}




