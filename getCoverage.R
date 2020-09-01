# Rui Shao
# 2020 Apr

options(scipen=999)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, GenomicRanges, IRanges, dplyr, rtracklayer,
               GenomicAlignments, BSgenome, Rsamtools, foreach, doParallel)

registerDoParallel(cores = 4)

processCoverage <- function(file_name, which_genome = 'mm10')
{ # this function converts input bam/bw file to coverage rle list
  # if input a paired-end bam file, coverage combines the reads strandness
  suffix <- toupper(sapply(strsplit(file_name,'\\.'), tail, 1))

  if (suffix == 'BW')
  {
    chr_lengths <- BSgenome::getBSgenome(which_genome) %>% seqlengths

    input_cov <- rtracklayer::import.bw(file_name)

    .seqnames <- as.character(seqnames(input_cov))
    .seqnames <<- .seqnames[ .seqnames %in% names(chr_lengths) ]

    input_cov <- foreach( chr = unique(.seqnames) ) %dopar%
    {
      tmp.gr = input_cov[ .seqnames == chr ]
      coverage( IRanges( start = start(tmp.gr), width = width(tmp.gr) ),
                weight = tmp.gr$score, width = chr_lengths[chr] )
    }
    names(input_cov) = unique(.seqnames)
  }
  if (suffix == 'BAM')
  {
    input_cov <- GenomicAlignments::readGAlignments(file_name) %>% coverage
  }
  return(input_cov)
}

binSums <- function(tmp_cov, bin_width, add_name = F)
{
  cov_len <- length(tmp_cov)
  i <- seq_len(cov_len %/% bin_width) * bin_width
  if (cov_len %% bin_width) i <- c(i, cov_len)
  bin_count <- diff( c(0L, cumsum(tmp_cov)[i]) )
  if (add_name) names(bin_count) <- i
  return(bin_count)
}

covMatrix <- function(intervals, input_cov, bin_width, new_len)
{ # Args:
  # intervals: genomic intervals to extract coverage matrix
  # input_cov: a list object from processCoverage()
  # extract coverage from input interval, save as a matrix
  .seqnames <- intersect.Vector(as.character(seqnames(intervals)), names(input_cov))

  out_matrix <- foreach( i = unique(.seqnames), .combine = cbind) %dopar%
  {
    tmp_gr <- intervals[.seqnames==i]
    tmp_cov <- input_cov[[i]]
    chr.end <- length(tmp_cov)

    sapply(seq_along(tmp_gr), function(x)
    {
      if(end(ranges(tmp_gr[x])) > chr.end)
      {
        tmp_binned <- c(tmp_cov, rep(0, end(ranges(tmp_gr[x])) - chr.end + 1))[ranges(tmp_gr[x])] %>%
          as.numeric %>% binSums(bin_width)
      } else {
        tmp_binned <- tmp_cov[ranges(tmp_gr[x])] %>% as.numeric %>% binSums(bin_width)
      }
      if (!is.na(new_len))
      {
        tmp_binned = smooth.spline(x = seq_along(tmp_binned), y = as.numeric(tmp_binned), df = 20)$y
        tmp_binned = spline(x=seq_along(tmp_binned), y=tmp_binned, n = new_len, method = 'natural')$y
        tmp_binned[tmp_binned < 0] = 0
      }
      tmp_binned
    })
  }
  out_matrix <- t(out_matrix)
  if (as.character(strand(intervals[1])) != '*')
  {
    neg_strand <- which( as.character(strand(intervals)) == '-' )
    out_matrix[neg_strand, ] <- t(apply(matrix(out_matrix[neg_strand, ],
                                               nrow = length(neg_strand)), 1, rev))
  }
  return(out_matrix)
}

## read file list
convertCoverage <- function(file_names, intervals, which_genome,
                            bin_width = 50, save_genome_cov = F, new_len = NA)
{ # Args:
  # file_names: BW file or BAM file without read strandness
  # intervals: GRanges intervals
  # save_genome_cov: if TRUE the full genome coverage Rle object will be in outputs
  # new_len: if NA this function will bin the genome
  # convert bw or bam files to a coverage array on input intervals

  if (is.na(new_len))
  {
    out_array <- array(dim = c(length(intervals),
                               ceiling(width(intervals[1])/bin_width),
                               length(file_names)))
  } else {
    out_array <- array(dim = c(length(intervals),
                               new_len,
                               length(file_names)))
  }

  gcov_list <- list()
  for ( i in seq_along(file_names))
  {
    input_cov <- processCoverage(file_names[i], which_genome)
    out_array[, , i] <- covMatrix(intervals, input_cov, bin_width, new_len)
    if(save_genome_cov) gcov_list <- c(gcov_list, list(input_cov))
  }

  if(save_genome_cov)
  {
    return(list(out_array = out_array, gcov_list = gcov_list))
  } else
  {
    return(out_array)
  }
}

binBwBam <- function(file_names, which_genome = 'mm10', bin_width = 1000L)
{
  # this function convert bw/bam files to binned matrix
  bin_dat <- NULL
  for (i in seq_along(file_names))
  {
    input_cov <- processCoverage(file_names[i], which_genome = which_genome)
    tmp_dat <- NULL
    for(j in names(input_cov))
    {
      tmp_bin_counts <- input_cov[[j]] %>% as.numeric %>% binSums(bin_width = bin_width, add_name = T)
      names(tmp_bin_counts) <- paste(j, names(tmp_bin_counts), sep = ": ")
      tmp_dat <- c(tmp_dat, tmp_bin_counts)
    }
    bin_dat <- cbind(bin_dat, tmp_dat)
  }
  colnames(bin_dat) <- gsub(".*\\/(.*)\\.*", "\\1", file_names)
  return(bin_dat)
}

readCoverage <- function(bam.files, object.gr, targets = "all", flank.size, is.frag = F)
{
  # count reads coverage with given bamfiles
  # output a list of Rle strings of each range
  # Agrs:
  #     bam.file: paired end strand specific bam file path
  #     object.gr: multiple intervals of gene or transcript
  #     targets: "full", "exon", "intron" (depends on Gene_input: complete gene reference)
  #     flank.size: flank both ends to find precise boundary or for plotting
  #     Gene_input: complete gene reference, GENCODE/ensembl
  #     is.frag: if coverage piles from paired-end reads fragments
  # refer to NGS.plot.r
  #
  object.gr = object.gr + flank.size
  object.start = start(object.gr)
  object.end = end(object.gr)
  object.strand = as.character(strand(object.gr)) == "+"

  if(targets != "all")
    exon.gr = Gene_input[Gene_input$type=="exon"]
    # exon.gr = Gene_input[Gene_input$type=="exon" & Gene_input$transcript_id%in%object.gr$transcript_id]

  for(i in seq_along(bam.files))
  {
    file.name = gsub('*.bam','\\',strsplit(bam.files[i],'\\/')%>%unlist%>%tail(1))
    file.dir = gsub(paste0('*', file.name,'.bam'), '\\', bam.files[i])
    bam.index = paste0(file.dir,file.name,'.bam.bai')
    if(!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.files[i])
    invisible(capture.output(is.paired.end <- testPairedEndBam(bam.files[i], index = bam.index)))

    .is.frag = is.frag
    if(!is.paired.end & is.frag) .is.frag = F
    # scanBamWhat: the info that need to be extracted from a bam file.
    sbw <- c('pos', 'qwidth', 'mapq', 'strand', 'rname',
             'mrnm', 'mpos', 'isize')
    sbp <- ScanBamParam(what = sbw, which = object.gr,
                        flag = scanBamFlag(isUnmappedQuery = F,
                                           isSecondaryAlignment = F,
                                           isFirstMateRead = is.paired.end))#use first in pair reads for coverage

    # Scan bam file to retrieve short reads.
    sr.in.ranges <- scanBam(bam.files[i], param = sbp, index = bam.index)

    combind.list = function(a, b)
    {
      new = lapply(names(a), function(x) c(a[[x]],b[[x]]))
      names(new) = names(a)
      new
    }

    scanBamRevOrder <- function(org.gr, sbp) {
      # ScanBamParam re-arranges the input genomic ranges. Use range info to
      # construct a string vector to find the order to reverse it.
      org.grnames <- with(org.gr, paste(seqnames, start, end, sep=':'))
      sbw.gr <- as.data.frame(bamWhich(sbp))  # scan-bam-ed
      if('space' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(space, start, end, sep=':'))
      } else if('group_name' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(group_name, start, end, sep=':'))
      } else {
        stop("Cannot locate chromosome names in extracted short reads. Report
             this problem using issue tracking or discussion forum.\n")
      }
      match(org.grnames, sbw.grnames)
    }

    # Restore the original order.
    sr.in.ranges <- sr.in.ranges[scanBamRevOrder(as.data.frame(object.gr), sbp)]

    srg.temp = list()
    for(n in seq_along(object.gr))
    {
      # Filter short reads by mapping quality and insertion size.
      all.index <- with(sr.in.ranges[[n]], (mapq>3 & strand == as.character(strand(object.gr[n])) ))
      if(is.paired.end)
      {
        all.index <- all.index & with(sr.in.ranges[[n]], !is.na(isize))
        sr.in.ranges[[n]]$isize = with(sr.in.ranges[[n]], mpos-pos )
      }
      srg.temp <- c(srg.temp, list(lapply(sr.in.ranges[[n]], function(x) x=x[all.index] )) )
    }
    names(srg.temp) = names(sr.in.ranges)
    # combine multiple bam files
    if(!exists('srg.filtered'))
    {
      srg.filtered = srg.temp
    }else
    {
      srg.filtered = lapply(names(srg.temp),
                            function(y) new = combind.list(srg.temp[[y]], srg.filtered[[y]]) )
    }
  }
  cov.list = list()
  for(n in seq_along(object.gr))
  {
    if(targets == "all")
    {
      cov = with( srg.filtered[[n]],
                  coverage(IRanges(start = ifelse(rep(object.strand[n], length(pos)), pos, mpos) - object.start[n],
                                   width = ifelse(rep(.is.frag, length(isize)), abs(isize), 0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    } else if(targets == "exon")
    {
      ov.index = with( srg.filtered[[n]],
                       findOverlaps(IRanges(start = pos, width = qwidth),
                                                       ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                         countQueryHits() != 0 )

      cov = with( lapply(srg.filtered[[n]], function(x) x[ov.index]),
                  coverage(IRanges(start = pos - object.start[n],
                                   width = ifelse(.is.frag, abs(isize),0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    } else if(targets == "intron")
    {
      ov.index = with( srg.filtered[[n]],
                       findOverlaps(IRanges(start = pos, width = qwidth),
                                    ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                         countQueryHits()==0 )

      cov = with( lapply(srg.filtered[[n]], function(x) x[ov.index]),
                  coverage(IRanges(start = pos - object.start[n],
                                   width = ifelse(.is.frag, abs(isize),0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    }

    if(!object.strand[n]) cov = rev(cov)
    cov.list = c(cov.list, cov)
  }
  return(cov.list)
}

resizeCov = function(cov.list, df=10, len=100)
{
  foreach(i=seq_along(cov.list), .combine = rbind) %dopar%
  {
    if (class(cov.list[i]) == "list") cov.tmp = cov.list[[i]] else cov.tmp = cov.list[i]
    cov = smooth.spline(x = seq_along(cov.tmp), y = as.numeric(cov.tmp), df = df)$y
    cov = spline(x = seq_along(cov), y = cov, n = len, method = 'natural')$y
    cov[cov < 0] = 0
    cov
  }
}

# bamsignals::bamProfile(bampath = bam.files, gr = object.gr, ss = T, paired.end = "ignore")

# read counting and coverage counting ------------------------------------------------------

.countBW <- function(bw_files, intervals, fast = F)
{ # Args:
  # extract interval coverage from unstranded bw files
  # save as a table,  (intervals, bw_files)

  if (fast) {
    count_dat <- sapply(bw_files, function (file)
      unlist(summary(BigWigFile(file), intervals))$score) %>%
      as.data.frame
  } else {
    count_dat <- foreach ( i = seq_along(bw_files), .combine = cbind) %dopar%
      {
        # bw_i <- rtracklayer::import(bw_files[i], which = intervals, format = 'BW')
        # matches <- findOverlaps(bw_i, intervals, type = "any")
        # aggregate(score(bw_i)[queryHits(matches)], by = list(subjectHits(matches)), sum)[, 2]
        bw_i <- rtracklayer::import(bw_files[i], format = 'BW', as="Rle")
        interval_counts <- NULL
        for(i in seq_along(intervals))
        {
          tmp.chr <- seqnames(intervals[i]) %>% as.character
          interval_counts <- c(interval_counts,
                               sum(Views(bw_i[[tmp.chr]], ranges(intervals[i]))) )
        }
        return(interval_counts)
      } %>% as.data.frame
  }
  colnames(count_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bw_files)
  rownames(count_dat) <- names(intervals)
  return(count_dat)
}

.countBam <- function(bam_files, intervals, stranded = F, paired.end = "ignore")
{ # extract count reads from bam files
  count_dat <- foreach ( bam_path = bam_files, .combine = cbind) %dopar%
  {
    if(!file.exists(paste0(bam_path,'.bai'))) indexBam(bam_path)
    read_counts <- bamsignals::bamCount(bam_path,
                                        intervals,
                                        paired.end = paired.end,
                                        ss = stranded, # split strandness
                                        verbose = F)
    if(stranded) read_counts <- read_counts[1, ]
    return(read_counts)
  } %>% as.data.frame
  colnames(count_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  if(is.null(names(intervals))) {
    rownames(count_dat) <- seq_along(intervals)
  } else {
    rownames(count_dat) <- names(intervals)
  }
  return(count_dat)
}

.coverBam <- function(bam_files, intervals, paired.end = "ignore", df = 20, len = 100)
{ # extract base-wise coverage of bam files reads
  # paired.end != "ignore", only first read in pairs
  # paired.end=="extend", fragment is treated as a single read
  registerDoParallel(cores = 2)
  cov_list <- foreach ( bam_path = bam_files) %dopar%
  {
    if(!file.exists(paste0(bam_path,'.bai'))) indexBam(bam_path)
    cov.list = bamsignals::bamCoverage(bam_path,
                                       gr = intervals,
                                       paired.end = paired.end,
                                       verbose = F)
    resizeCov(cov.list, df = df, len = len)
  }
  names(cov_list) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  return(cov_list)
}

extend_range <- function(x, up_flank, down_flank)
{
  .strand = as.character(strand(x))
  start(x) = start(x) - ifelse(.strand %in% c("+", "*"), up_flank, down_flank)
  end(x) = end(x) + ifelse(.strand %in% c("-", "*"), up_flank, down_flank)
  x
}

readBam <- function(bam_files, intervals, pair_end = F, stranded = F,
                    flanks = NULL, # e.g. flanks = c(up=2000, down=2000),
                    new_lens = c(up=100, mid=50, down=200))
{ # extract interval read coverage with flanking regions from bam files
  # save to a list object for each bam file
  if (is.null(flanks) & length(unique(width(intervals))) > 1) stop("Interval without resizing must be the same width.")
  registerDoParallel(cores = 10)
  
  intervals <- keepSeqlevels(intervals, as.character(unique(seqnames(intervals))))
  
  if (!is.null(flanks)) 
  {
    intervals <- extend_range(intervals, flanks[1], flanks[2])
  } else {
    flanks <- c(0, 0)
  }
  
  cov_dat <- foreach(bam_path = bam_files) %dopar%
  {
      bam.index=paste0(bam_path,'.bai')
      if(!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam_path)
      if(pair_end)
      {
        sbw = c('pos', 'qwidth','strand','rname', 'mrnm', 'mpos', 'isize')
        flag = scanBamFlag(isFirstMateRead = pair_end, isSecondaryAlignment=F)
      } else {
        sbw = c('pos', 'qwidth','strand','rname')
        flag=scanBamFlag(isSecondaryAlignment=F)
      }
      param = ScanBamParam(what = sbw, flag = flag, which = intervals)
      srg = scanBam(bam_path, param=param, index = bam.index)
      .cov_mat = NULL
      for(n in seq_along(intervals))
      {
        .strand = as.character(strand(intervals[n]))
        if(stranded)
        {
          all.index <- with(srg[[n]],
                            !is.na(mpos) & mrnm == rname & abs(isize) < 2000 & strand == .strand)
        } else {
          all.index <- with(srg[[n]], !is.na(mpos) & mrnm == rname & abs(isize) < 2000 )
        }

        srg.tmp <- lapply(srg[[n]], function(x) x=x[all.index])
        .cov = with(srg.tmp, coverage(IRanges(start = pos - start(intervals[n]),
                                              width = qwidth ),
                                      width = width(intervals[n]) ))
        if ( .strand == '-') .cov = rev(.cov)
        .cov = smooth.spline(x = seq_len(width(intervals[n])), y=as.numeric(.cov), df = 200)$y

        if (sum(flanks) > 0)
        {
          .cov_resize = NULL
          if (flanks[1] > 0) .cov_resize = spline(x=seq_len(flanks[1]),
                                                  y=.cov[seq_len(flanks[1])],
                                                  n = new_lens[1],
                                                  method = 'natural')$y
          
          .cov_resize = c(.cov_resize, spline(x = seq_len(width(intervals[n]) - sum(flanks)),
                                              y = .cov[(sum(flanks[1]) + 1) : (width(intervals[n]) - sum(flanks[2]))],
                                              n = new_lens[2],
                                              method = 'natural')$y)
          
          if (flanks[2] > 0) .cov_resize = c(.cov_resize, spline(x=seq_len(sum(flanks[2])),
                                                                      y=.cov[(width(intervals[n]) - sum(flanks[2]) + 1) : width(intervals[n])],
                                                                      n = new_lens[3],
                                                                      method = 'natural')$y)
        }
        .cov_resize[.cov_resize < 0] = 0
        .cov_mat = rbind(.cov_mat, .cov_resize)
      }
      rownames(.cov_mat) <- names(intervals)
      .cov_mat
    } # end of processing all bam files
  names(cov_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  return(cov_dat)
}

intervalCoverage <- function(query.list, intervals, out_width)
{ # this function converts GRange objects in query.list to coverage matrix upon subject intervals
  # e.g. find motif coverage on intervals
  cov.list <- list()

  for(query in query.list)
  {
    cov.mat.p <- matrix(0, nrow = length(intervals), ncol = out_width)
    cov.mat.m <- matrix(0, nrow = length(intervals), ncol = out_width)

    mtch.p <- findOverlaps(query, intervals, ignore.strand = F) # plus overlap
    levels(strand(query)) <- c('-', '+', '*') %>% as.factor
    mtch.m <- findOverlaps(query, intervals, ignore.strand = F) # minus overlap

    s.p <- split(queryHits(mtch.p), list(subjectHits(mtch.p)))
    s.m <- split(queryHits(mtch.m), list(subjectHits(mtch.m)))

    for (i in names(s.p))
    {
      i.n <- as.numeric(i)
      start.p <- start(query[ s.p[[i]] ]) - start(intervals[i.n])
      width.p <- width(query[ s.p[[i]] ])

      cov.mat.p[i.n, ] <- coverage( IRanges( start = start.p, width = width.p ),
                                    width = width(intervals[i.n]) ) %>% spline(n=out_width) %>% '$'(y)
    }

    for ( i in names(s.m))
    {
      i.n <- as.numeric(i)
      start.m <- start(query[ s.m[[i]] ]) - start(intervals[i.n])
      width.m <- width(query[ s.m[[i]] ])

      cov.mat.m[i.n, ] <- coverage( IRanges( start = start.m, width = width.m ),
                                    width = width(intervals[i.n]) ) %>% spline(n=out_width) %>% '$'(y)
    }

    minus.intevals <- as.character(strand(intervals)) == '-'
    cov.mat.p[minus.intevals, ] <- rev(cov.mat.p[minus.intevals, ])
    cov.mat.m[minus.intevals, ] <- rev(cov.mat.m[minus.intevals, ])

    cov.list <- c(cov.list, list(list('plus' = cov.mat.p, 'minus' = cov.mat.m)))
  }
  names(cov.list) <- (names(query.list))
  return(cov.list)
}

viewCoverage <- function(bam.file.iist, bam.sample.names = NULL,
                         interval, 
                         stranded = T, is.fragment = F,
                         size.factor.list = NULL, 
                         log_scale = T, low_cut = 2,
                         smoothen = F, df = 300,
                         bin_width = 1,
                         cov_color_list = NULL,
                         annotation_list = NULL, 
                         anno_box_height = 2, 
                         anno_text_cex = 1.5) {
  
  # Initialize parameters
  if (is.null(names(bam.file.iist))) 
    names(bam.file.iist) <- seq_along(bam.file.iist)
  if (is.null(bam.sample.names)) 
    bam.sample.names <- names(bam.file.iist)
  if (is.null(size.factor.list)) 
    size.factor.list <- lapply(bam.file.iist, function(x) rep(1, length(x)))
  if (any(names(size.factor.list) != names(bam.file.iist))) 
    names(size.factor.list) <- names(bam.file.iist)
  if (is.null(cov_color_list)) 
    cov_color_list <- lapply(bam.file.iist, function(x) c("plus"='#bd5734', "minus"='#feb236'))
  
  combine_rep <- function(cov.list, size.factors, log_scale, 
                          smoothen = F, df = 300) { 
    # get the mean normalised coverage
    out <- Reduce("+", sapply(seq_along(cov.list), 
                              function(i) cov.list[[i]] / size.factors[i]))
    out <- as.numeric(out / length(cov.list))
    if (smoothen) out <- smooth.spline(seq_along(out), y = out, df = df)$y
    out[out < 0] <- 0
    if (log_scale) {
      out[out < low_cut] <- 1
      return(log10(out))
    }
    out
  }
  
  # Initialize an empty list
  coverage_list <- vector(mode = "list", length = length(bam.file.iist))
  names(coverage_list) <- names(bam.file.iist)
  
  # Get coverage and combine replicates
  for (sample in names(bam.file.iist)) 
  {
    bam_files = bam.file.iist[[sample]]
    plus.cov = 
      combine_rep(
        cov.list = unlist(lapply(bam_files, readCoverage, object.gr = interval, flank.size=0, is.frag = is.fragment)),
        size.factors = size.factor.list[[sample]], 
        log_scale = log_scale, 
        smoothen = smoothen, 
        df = df)
    if (stranded) 
    {
      levels(strand(interval)) = c("-", "+", "*")
      minus.cov = 
        combine_rep(
          cov.list = unlist(lapply(bam_files, readCoverage, object.gr = interval, flank.size=0, is.frag = is.fragment)),
          size.factors = size.factor.list[[sample]], 
          log_scale = log_scale, 
          smoothen = smoothen, 
          df = df)
      minus.cov = rev(minus.cov)
      levels(strand(interval)) = c("-", "+", "*")
      # save interval coverage into list for plotting
      coverage_list[[sample]] = list("plus" = binSums(plus.cov, bin_width = bin_width) / bin_width, 
                                     "minus" = binSums(minus.cov, bin_width = bin_width) / bin_width)
    } else {
      coverage_list[[sample]] = list("plus" = binSums(plus.cov, bin_width = bin_width) / bin_width)
    }
  }
  
  # Plot start
  par(mfrow = c(length(coverage_list) + length(annotation_list), 1))
  par(mar = c(0,5,1,1))
  
  x_lims = c(0, length(coverage_list[[1]][[1]]))
  if (stranded) {
    max(unlist(lapply(coverage_list, function(x) max(x[["minus"]]))))
    y_lims = c(-max(unlist(lapply(coverage_list, function(x) max(x[["minus"]])))),
               max(unlist(lapply(coverage_list, function(x) max(x[["plus"]])))))
  } else {
    y_lims = c(0, max(unlist(lapply(coverage_list, function(x) max(x[["plus"]])))))
  }
  y_off_set = y_lims[2] / 30
  
  for (sample in names(bam.file.iist)) 
  {
    plot(x_lims, y_lims, type = "n",bty="n", xaxt='n', xlab="", yaxt = "n", ylab = sample, cex.lab = 1.2)
    axis(side = 2, lwd=1.5, cex.axis=1.2)
    barplot(height = coverage_list[[sample]][["plus"]], col=cov_color_list[[sample]][1],
            add = TRUE, axes = FALSE, space = 0, border = NA, offset = y_off_set)
    if (stranded)
      barplot(height = -coverage_list[[sample]][["minus"]], col=cov_color_list[[sample]][2], 
              add = TRUE, axes = FALSE, space = 0, border = NA, offset = -y_off_set)
  }
  
  for (anno_name in names(annotation_list)) {
    # Extract annotations in the plotting window
    if (length(annotation_list[[anno_name]] == 0)) next()
    anno = annotation_list[[anno_name]]
    anno = anno[countQueryHits(findOverlaps(anno, interval, ignore.strand = T)) > 0]
    
    if (anno_name == "Gene") {
      exons = reduce(anno[anno$type == "exon"]) 
      genes = anno[anno$type == "gene"]
      genes.plus.anno = genes[strand(genes) == "+"]
      genes.minus.anno = genes[strand(genes) == "-"]
      
      exon.plus = intersect(IRanges(start = 1, width = width(interval)), 
                            shift(ranges(exons[as.character(strand(exons)) == "+"]),
                                  shift = -start(interval)))
      exon.minus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(exons[as.character(strand(exons)) == "-"]),
                                   shift = -start(interval)))
      genes.plus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(genes[as.character(strand(genes)) == "+"]),
                                   shift = -start(interval)))
      genes.minus = intersect(IRanges(start = 1, width = width(interval)), 
                              shift(ranges(genes[as.character(strand(genes)) == "-"]),
                                    shift = -start(interval)))
      # Height for plotting
      h.exon.plus = h.exon.minus = h.intron.plus = h.intron.minus = Rle(0, width(interval))
      h.exon.plus[exon.plus] = anno_box_height
      h.exon.minus[exon.minus] = anno_box_height
      h.intron.plus[genes.plus] = anno_box_height / 4
      h.intron.minus[genes.minus] = anno_box_height / 4
      
      # gene ranges
      start.plus = start(genes.plus)
      start.minus = end(genes.minus)
      mid.plus = start(resize(genes.plus, fix = "center", 1))
      mid.minus = start(resize(genes.minus, fix = "center", 1))
      
      # Gene boxes
      barplot(height = as.numeric(h.intron.plus), col='#0000B2', ylim = c(-10, 10), add = FALSE,
              axes = FALSE, space = 0, border = NA, offset = 3.7)
      barplot(height = -as.numeric(h.intron.minus), col='#0000B2', add = TRUE,
              axes = FALSE, space = 0, border = NA, offset = -3.7)
      barplot(height = as.numeric(h.exon.plus), col='#0000B2', add = TRUE,
              axes = FALSE, space = 0, border = NA, offset = 3)
      barplot(height = -as.numeric(h.exon.minus), col='#0000B2', add = TRUE,
              axes = FALSE, space = 0, border = NA, offset = -3)
      
      
      # height <- 1
      # if (is(xlim, "IntegerRanges"))
      #   xlim <- c(min(start(xlim)), max(end(xlim)))
      # bins <- disjointBins(IRanges(start(x), end(x) + 1))
      # plot.new()
      # plot.window(xlim, c(0, max(bins)*(height + sep)))
      # ybottom <- bins * (sep + height) - height
      # rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
      
      
      # gene names
      label.plus = genes.plus.anno$gene_name
      label.minus = genes.minus.anno$gene_name
      
      # remove embeded genes
      if (length(genes.plus) < sum(strand(genes) == "+")) 
        label.plus = genes.plus.anno$gene_name[countQueryHits(findOverlaps(genes.plus.anno, genes.plus.anno, type = "within")) < 2]
      if (length(genes.minus) < sum(strand(genes) == "-")) 
        label.minus = genes.plus.anno$gene_name[countQueryHits(findOverlaps(genes.minus.anno, genes.minus.anno, type = "within")) < 2]
      
    } else if (anno_name == "TU") {
      
      anno.plus = anno[strand(anno) == "+"]
      anno.minus = anno[strand(anno) == "-"]
      
      range.plus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(anno.plus), shift = -start(interval)))
      range.minus = intersect(IRanges(start = 1, width = width(interval)), 
                              shift(ranges(anno.minus), shift = -start(interval)))
      
      if (length(range.plus) < length(anno.plus)) anno.plus = anno.plus[anno.plus$location %ni% c("dsRNA", "usRNA")]
      if (length(range.minus) < length(anno.minus)) anno.minus = anno.minus[anno.minus$location %ni% c("dsRNA", "usRNA")]
      
      # Height for plotting
      cov.plus = cov.minus = Rle(0, width(interval))
      cov.plus[range.plus] = anno_box_height
      cov.minus[range.minus] = anno_box_height
      
      # TU ranges
      start.plus = start(range.plus)
      start.minus = end(range.minus)
      mid.plus = start(resize(range.plus, fix = "center", 1))
      mid.minus = start(resize(range.minus, fix = "center", 1))
      
      # TU rectangular
      barplot(height = as.numeric(cov.plus), col='#bc1066', ylim = c(-10, 10), add = FALSE,
              axes = FALSE, space = 0, border = NA, offset = 3)
      barplot(height = -as.numeric(cov.minus), col='#bc4066', add = TRUE,
              axes = FALSE, space = 0, border = NA, offset = -3)
      
      # TU names
      label.plus = anno.plus$location
      label.minus = anno.minus$location
    }
    
    # TSS arrow and TU name
    if (length(label.plus) > 0) 
    {
      text(x = mid.plus, y = 7, labels = label.plus, cex = anno_text_cex)
      arrows(x0 = start.plus, y0 = 8, 
             x1 = start.plus + width(interval) / 50, y1 = 8, 
             code = 2, angle = 10, cex = 0.1, length = 0.05)
      arrows(x0 = start.plus, y0 = 3.5, 
             x1 = start.plus, y1 = 7.9, length = 0)
    }
    
    if (length(label.minus) > 0) 
    {
      text(x = mid.minus, y = -7, labels = label.minus, cex = anno_text_cex)
      arrows(x0 = start.minus, y0 = -8, 
             x1 = start.minus - width(interval) / 50, y1 = -8, 
             code = 2, angle = 10, cex = 0.1, length = 0.05)
      arrows(x0 = start.minus, y0 = -3.5, 
             x1 = start.minus, y1 = -7.9, length = 0)
    }
  }
}
