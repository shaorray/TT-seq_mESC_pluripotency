# utilities

pacman::p_load(BiocManager, rtracklayer,Rsamtools,
               GenomicRanges, IRanges, org.Mm.eg.db,
               foreach, doParallel, 
               dplyr, tidyverse, tibble, matrixStats,
               ggplot2, grid, gridExtra, viridis, RColorBrewer)

registerDoParallel(cores = 4)

importRanges <- function(data, seq_lvls = paste0('chr',c(1:21,'X','Y')))
{
  format <- toupper(sapply(strsplit(data,'\\.'), tail, 1))
  import_input <- rtracklayer::import(data, format = format)
  seqlevelsStyle(import_input) <- "UCSC"
  import_input <- import_input[ seqnames(import_input) %in% seq_lvls ]
  seqlevels(import_input) <- as.character(unique(seqnames(import_input)))
  return(import_input)
}

readKallistoResults <- function( file_paths, sample_names )
{
  txCounts <- NULL
  for(i in file_paths)
  {
    temp.table <- read.table(file = i, header=T, colClasses=c("character", "numeric", "numeric", "numeric", "numeric"))
    txCounts <- cbind(txCounts, temp.table[, 4])
  }
  rownames(txCounts) <- unname(temp.table$target_id)
  colnames(txCounts) <- sample_names[seq_along(file_paths)]
  txCounts <- keepOneTx(txCounts) # the max in tx variants
  return(txCounts)
}

keepOneTx <- function(count_table, rowname_gene_id = F)
{ # process kallisto counts on genecode fasta reference
  # save only the main tx counts of protein_coding and lincRNA genes

  trim_type <- function(x) gsub("(.*\\|)(.*)(\\|$)", "\\2", x)
  gene_types <- rownames(count_table) %>% trim_type
  count_table <- count_table[gene_types %in% c('protein_coding', 'lincRNA'), ]
  out_tx_id <- T
  trim_gene <- function(x) lapply(strsplit(x, '\\|'),
                                  function(y) gsub("*(\\..*)",
                                                   "\\2",
                                                   y[ifelse(out_tx_id, 1, 2)])
                                  ) %>% unlist
  tx_ids <- rownames(count_table) %>% trim_gene
  out_tx_id <- F
  gene_ids <- rownames(count_table) %>% trim_gene
  idx <- split(tx_ids, list(gene_ids))

  rownames(count_table) <- tx_ids

  row_sums <- rowSums(count_table)
  doParallel::registerDoParallel(cores = 8)

  new_table <- foreach(i = idx, .combine = rbind) %dopar%
    {
      if(length(i) == 1)
      {
        matrix(count_table[i, ],
               nrow = 1,
               dimnames = list(i, colnames(count_table)))
      } else {
        tmp <- data.frame(count_table[i, ])
        as.matrix(tmp[which.max(row_sums[i]), ])
      }
    }

  colnames(new_table) <- colnames(count_table)
  if (rowname_gene_id) rownames(new_table) <- names(idx)

  return(new_table)
}

SizeFactorCal <- function(gene.counts, RefExpr=NULL)
{
  #estimate relative size factor comparing to reference TU gene counts
  #Args:
  #gene.counts: genes by row and samples by column
  #RefExpr: if defined the gene.counts table will return relative size factors
  gene.counts = gene.counts[!apply(gene.counts,1,function(x) any(is.na(x)) ) & rowSums(gene.counts)>0,]
  if(!is.null(RefExpr)){
    gene.counts = cbind(RefExpr, gene.counts)
    # geometric mean
    geoMeans = apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm=T) / length(x)) )
    # quotients list to median (size factor)
    quoMedian = apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median)
    return(quoMedian[-1] / quoMedian[1])
  }else{
    geoMeans=apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm=T) / length(x)) )
    return(apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median))
  }
}

rowMax <- function(mat) apply(mat, 1, max)

geoMeans <- function(mat)
{
  if (is.null( dim(mat) )) {
    mat
  } else {
    exp( log(rowProds(mat)) / ncol(mat))
  }
}

groupConsec <- function( seq ) # -> order
{ # input a sequence of classes
  # sum up consecutive numbers
  vals <- seq[1]
  val_group <- 1
  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  for ( i in seq_along(seq)[-1] )
  {
    if (seq[i-1] == seq[i])
    {
      val_group[length(val_group)] %+=% 1
    } else
    {
      vals <- c( vals, seq[i] )
      val_group <- c(val_group, 1)
    }
  }
  names(val_group) <- vals
  return(val_group)
}

kink_index <- function(x, method = "slope")
{
  # this function returns outlier indexes of the input vector
  # to retrieve highly enriched loci, e.g. HOMER uses slope = 1 to get super enhancers
  x[is.na(x)] <- 0
  x_ord <- order(x)
  if (method == "slope")
  {
    slopes_above_step <- diff(x[x_ord]) > 1
    stopifnot(sum(slopes_above_step) > 0)
    return(x_ord[seq(min(which(slopes_above_step)) + 1, length(x))])
  }
  if (method == "weight")
  {
    x[is.infinite(x)] <- max(x[!is.infinite(x)])
    weight_above_mean <- which.min(cumsum(diff(x[x_ord]) - mean(diff(x[x_ord]))))
    return(x_ord[seq(weight_above_mean, length(x))])
  }
}

plot_Vennerable <- function (list_1, list_2, 
                             name_1, name_2,
                             color_set = c(1, 2, 3), 
                             color_Text = c(1, 2), ...)
{
  Sets = list(list_1, list_2)
  names(Sets) = c(name_1, name_2)
  p = Vennerable::compute.Venn(Vennerable::Venn(Sets = Sets))
  gp = Vennerable::VennThemes(p)
  gp$Face$`10`$fill = color_set[1]
  gp$Face$`11`$fill = color_set[2]
  gp$Face$`01`$fill = color_set[3]
  gp$FaceText$`10`$fontsize = 30
  gp$FaceText$`11`$fontsize = 30
  gp$FaceText$`01`$fontsize = 30
  gp$Set$Set1$lwd = 0
  gp$Set$Set2$lwd = 0
  gp$SetText$Set1$fontsize = 32
  gp$SetText$Set2$fontsize = 32
  gp$SetText$Set1$col = color_Text[1]
  gp$SetText$Set2$col = color_Text[2]
  # plot and add pval
  plot(p, gp = gp, show = list(Universe=F), ...)
  pacman::p_load(GeneOverlap)
  grid.text(
    paste(
      # "p <", formatC(testGeneOverlap(newGeneOverlap(list_1, list_2))@pval,
      #                format = "e", digits = 2),"; ",
      "Jaccard =", round(length(intersect.Vector(list_1, list_2)) / 
                           length(unique(c(list_1, list_2))), 2)
  ),
  x = 0.5, y=0.95, gp=gpar(cex=2))
}

get_dens <- function(X, Y, n.grid = 100) {
  # from https://slowkow.com/notes/ggplot2-color-by-density/
  dens <- MASS::kde2d(X, Y, n = n.grid)
  ix <- findInterval(X, dens$x)
  iy <- findInterval(Y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  # X_cut <- cut(X, seq(min(X), max(X), diff(range(X)) / n.grid ) )
  # Y_cut <- cut(Y, seq(min(Y), max(Y), diff(range(Y)) / n.grid ) )
  # den <- as.numeric(sqrt(table(X_cut)[X_cut] * table(Y_cut)[Y_cut]))
}

# graphic settings ---------------------------------------------------------------------------
theme_setting <- list(theme_minimal(),
                      theme(legend.position = 'right',
                            legend.justification = c(0, 1),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            # panel.border = element_rect(colour = "black", fill=NA, size=1),
                            plot.title = element_text(size = 14, face="bold", vjust = 0.5),
                            axis.text=element_text(size=12, face = "bold"),
                            axis.title=element_text(size=15,face="plain"),
                            legend.title = element_text(size=14, face="bold"),
                            legend.text = element_text(size=12, face="plain"),
                            axis.ticks.x = element_line(size = 1), 
                            axis.ticks.y = element_line(size = 1),
                            # axis.ticks.length = unit(20, "pt"),
                            panel.border = element_rect(colour = "black", fill=NA, size=1.2))
                      )

colors_9 <- c('#59C7EB', '#FEA090', '#9AA0A7', '#077187', '#0A9086', '#3E5496', '#E0607E', '#8E2043', '#EFDC60')

colors_20 <- c('#bf4040', '#cc8800', '#808000', '#408000', '#006622', '#2d8659', '#5c8a8a',
               '#0073e6', '#4d4dff', '#5353c6', '#4d0099', '#660080', '#602060', '#c6538c',
               '#99003d', '#cc5200', '#b32400', '#663300', '#b3b300', '#4d9900')

colors_n <- c('#c90000', '#c94600', '#c99700', '#6f9c00', '#009c56', '#00838f',
              '#0550b3', '#3e0080', '#560659', '#adadad')

add.alpha <- function(col, alpha=1)
{
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

`%+=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 + e2))
`%-=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 - e2))
`%*=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 * e2))
`%/=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 / e2))
`%c=%` = function(e1, e2) eval.parent(substitute(e1 <- c(e1, e2)))
`%ni%` = Negate(`%in%`)
