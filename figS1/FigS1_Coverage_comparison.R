# Compare nascent RNA-seq methods by exon coverage
# Rui Shao
# Aug 2019
 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")
registerDoParallel(cores = 5)

L.sample.list <- list(TTseq1 = list.files("/mnt/0E471D453D8EE463/TT_seq_data/bam_sl_2i_mm10/", "SL_.*bam$", full.names = T),
                      TTseq2 = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2018_Jan_TTseq_TX1072_WT", "0h.*bam$", full.names = T),
                      PROseq1 = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2016_Jesse_PROSeq", "bam$", full.names = T),
                      PROseq2 = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2018_Marta_ESC_PRO-Seq/", "bam$", full.names = T),
                      GROseq1 = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2016_Flynn_GRO-seq-mES-12hrKD_Control/", "bam$", full.names = T),
                      GROseq2 = list.files("/mnt/0E471D453D8EE463/GEO_nascent_RNA/2017_Dorighi_Gro-seq_WT/", "bam$", full.names = T))

Gene_input <- importRanges("/mnt/0E471D453D8EE463/genomeDir/Mus_musculus.GRCm38.79.gtf")
TU_input <- importRanges("../data/TU_anno/other/TU_filter+LRNA_SL_rep1.Aligned.sortedByCoord.out.gtf")

# get ExonIntron coverage
EI_res <- extractExonIntron(Gene_input, TU_input$transcript_id)

resSize = 600
cov.TEI.list <- coverageExonIntron(L.sample.list, EI_res, resSize = resSize, size.factor = size.factor)

# saveRDS(cov.TEI.list, "data/cov.TEI.list.RData")

colors=sapply(c(RColorBrewer::brewer.pal(10,'Paired'), RColorBrewer::brewer.pal(8,'Dark2')[6:7]),
              function(x) rgb(col2rgb(x)[1]/255,col2rgb(x)[2]/255,col2rgb(x)[3]/255,0.7) )

if (T)
{
  pdf("figs/FigS1_Exon_coverage.pdf", width = 6, height = 6)
  par(mfrow=c(3,2))
  par(mar=c(3,2,2,2))
  ###TT-seq###
  plot(c(0,resSize*3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  axis(side=2,lwd=1.5,cex.axis=1.2)
  barplot(height = cov.TEI.list$TSS[,3][seq_len(resSize * 3)],
          col=colors[1],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TSS[,3][seq_len(resSize * 3) + resSize * 3],
          col=colors[2],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = -cov.TEI.list$TSS[,4][seq_len(resSize * 3)],
          col=colors[2],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TSS[,4][seq_len(resSize * 3) + resSize * 3],
          col=colors[1],add=T,axes=F,space=0,border=NA,offset=0.08)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  plot(c(0,resSize*3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  barplot(height = cov.TEI.list$TTS[,3][seq_len(resSize * 3)],
          col=colors[1],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TTS[,3][seq_len(resSize * 3) + resSize * 3],
          col=colors[2],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = -cov.TEI.list$TTS[,4][seq_len(resSize * 3)],
          col=colors[2],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TTS[,4][seq_len(resSize * 3) + resSize * 3],
          col=colors[1],add=T,axes=F,space=0,border=NA,offset=0.08)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  ###PRO-seq###
  plot(c(0,resSize * 3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  axis(side=2,lwd=1.5,cex.axis=1.2)
  barplot(height = -cov.TEI.list$TSS[,5][seq_len(resSize * 3)],
          col=colors[4],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TSS[,5][seq_len(resSize * 3) + resSize * 3],
          col=colors[3],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TSS[,6][seq_len(resSize * 3)],
          col=colors[4],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TSS[,6][seq_len(resSize * 3) + resSize * 3],
          col=colors[3],add=T,axes=F,space=0,border=NA,offset=0.08)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  plot(c(0,resSize * 3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  barplot(height = -cov.TEI.list$TTS[,5][seq_len(resSize * 3)],
          col=colors[4],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TTS[,5][seq_len(resSize * 3) + resSize * 3],
          col=colors[3],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TTS[,6][seq_len(resSize * 3)],
          col=colors[4],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TTS[,6][seq_len(resSize * 3) + resSize * 3],
          col=colors[3],add=T,axes=F,space=0,border=NA,offset=0.08)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  ###GRO-seq###
  plot(c(0,resSize * 3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  axis(side=2,lwd=1.5,cex.axis=1.2)
  barplot(height = cov.TEI.list$TSS[,7][seq_len(resSize * 3)],
          col=colors[5],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TSS[,7][seq_len(resSize * 3) + resSize * 3],
          col=colors[6],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TSS[,8][seq_len(resSize * 3)],
          col=colors[5],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TSS[,8][seq_len(resSize * 3) + resSize * 3],
          col=colors[6],add=T,axes=F,space=0,border=NA,offset=-0.06)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  plot(c(0,resSize * 3),c(-3,5),type='n',bty="n",xaxt='n',xlab='',yaxt="n")
  barplot(height = cov.TEI.list$TTS[,7][seq_len(resSize * 3)],
          col=colors[5],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TTS[,7][seq_len(resSize * 3) + resSize * 3],
          col=colors[6],add=T,axes=F,space=0,border=NA,offset=-0.06)
  barplot(height = cov.TEI.list$TTS[,8][seq_len(resSize * 3)],
          col=colors[5],add=T,axes=F,space=0,border=NA,offset=0.08)
  barplot(height = -cov.TEI.list$TTS[,8][seq_len(resSize * 3) + resSize * 3],
          col=colors[6],add=T,axes=F,space=0,border=NA,offset=-0.06)
  abline(v=resSize,lwd=0.2,lty=2);abline(v=resSize*2,lwd=0.2,lty=2)
  
  dev.off()
}

# functions

extractExonIntron <- function(Gene_input, tx_ids)
{
  # ranges of first exon + first intron
  # and last intron + last exon
  # then 1000bp before TSS and after 1000bp TTS on both strand
  
  # TU annotation tx_ids restrict only active gene exons
  Gene_exon = Gene_input[Gene_input$type=='exon' & Gene_input$transcript_id %in% tx_ids]
  ExonIntron_1 = ExonIntron_2 = GRanges()
  range.index.list = list()
  exon.list = split(Gene_exon, Gene_exon$transcript_id)
  
  for (i in names(exon.list))
  { # select the first and the last exons
    exon.temp = exon.list[[i]] %>% sort
    E.1 = E.2 = exon.temp[1]
    
    if (length(exon.temp) >= 2)
    {
      if (as.character(strand(exon.temp[1])) == '+')
      {
        E.2 = exon.temp[length(exon.temp)] # last exon
        start(E.1) = start(E.1)-1000; end(E.1) = start(exon.temp[2]) # the end of 1st intron
        start(E.2) = end(exon.temp[length(exon.temp)-1]); end(E.2)=end(E.2) + 1000
        current.index.list = c(list(c(1,1000),
                                    c(1001, width(exon.temp[1])+1000),
                                    c(1001+width(exon.temp[1]), width(E.1) )),
                               list(c(1, start(exon.temp[length(exon.temp)])-end(exon.temp[length(exon.temp)-1])),
                                    c(start(exon.temp[length(exon.temp)])-end(exon.temp[length(exon.temp)-1])+1, width(E.2)-1000),
                                    c(width(E.2)-999, width(E.2) )))
      } else {
        E.1 = exon.temp[length(exon.temp)]
        end(E.1) = end(E.1)+1000; start(E.1) = end(exon.temp[length(exon.temp)-1])
        end(E.2) = start(exon.temp[2]); start(E.2) = start(exon.temp[1])-1000
        current.index.list = c(list(c(1,1000),
                                    c(1001, width(exon.temp[length(exon.temp)])+1000),
                                    c(1001+width(exon.temp[length(exon.temp)]), width(E.1)) ),
                               list(c(1, start(exon.temp[2])-end(exon.temp[1])),
                                    c(start(exon.temp[2])-end(exon.temp[1])+1, width(E.2)-1000),
                                    c(width(E.2)-999, width(E.2)) ) )
      }
      
      ExonIntron_1 = c(ExonIntron_1, E.1); ExonIntron_2 = c(ExonIntron_2, E.2)
      names(current.index.list) = c("TSS","Exon1","Intron1","Intron2","Exon2","TTS")
      current.index.list = list(current.index.list); names(current.index.list)=i
      range.index.list = c(range.index.list, current.index.list)
    } 
  }
  return(list(ExonIntron_1 = ExonIntron_1, ExonIntron_2 = ExonIntron_2,
              range.index.list = range.index.list))
}

coverageExonIntron <- function(L.sample.list, EI_res, resSize, size.factor = NULL)
{
  # get coverage of first exon + first intron
  # and last intron + last exon
  # then 1000bp before TSS and after 1000bp TTS on both strand
  
  # utils
  # .covResize = function(cov, len) spline(x=seq_along(cov), y=as.numeric(cov),n = len, method = 'natural')$y
  .covResize = function(cov, resSize) cov[ceiling(seq_len(resSize) * length(cov) / resSize)]
  # .covNorm = function(cov.mat) apply(cov.mat, 1, function(x) log(x+1)) %>% t
  .covSmoothen = function(cov, df=10) smooth.spline(x = seq_along(cov), y=cov, df = df)$y
  
  ExonIntron_1 = EI_res$ExonIntron_1
  ExonIntron_2 = EI_res$ExonIntron_2
  range.index.list = EI_res$range.index.list
  
  if (is.null(size.factor)) 
  {
    lapply(L.sample.list, function(x) 
      rowMeans(.countBam(bam_files = x, intervals = ExonIntron_1, stranded = T))) %>%
      Reduce(cbind, .) %>% SizeFactorCal -> size.factor
  }

  # initialize result table, TSS/Exon/Intron
  cov.TEI.list = list()
  cov.TEI.list$TSS = data.frame(position = rep(seq_len(resSize), 6),
                                strand = c(rep(1, resSize*3), rep(0, resSize*3)))
  cov.TEI.list$TTS = data.frame(position = rep(seq_len(resSize), 6),
                                strand = c(rep(1, resSize*3), rep(0, resSize*3)))
  # L.sample.list=first.sample.list
  
  for(i in seq_along(L.sample.list))
  {
    bam.files = L.sample.list[[i]]
    sizeFactor = size.factor[i]
    # sameStrand = !all.strandFlipped[i]
    sameStrand = T
    tempCov.sense = readCoverage(bam.files, object.gr = c(ExonIntron_1, ExonIntron_2), flank.size=0)
    tempCov.1 = tempCov.sense[seq_along(ExonIntron_1)]
    tempCov.2 = tempCov.sense[length(ExonIntron_2) + seq_along(ExonIntron_2)]
    rm(tempCov.sense)
    levels(strand(ExonIntron_1)) = levels(strand(ExonIntron_2)) = c("-", "+", "*") 
    cat('Sample+ '); cat(i); cat('\n')
    gc()
    tempCov.antisense = readCoverage(bam.files, object.gr = c(ExonIntron_1, ExonIntron_2),
                                     flank.size=0) %>% lapply(rev)
    tempCov.1a = tempCov.antisense[seq_along(ExonIntron_1)]
    tempCov.2a = tempCov.antisense[length(ExonIntron_2) + seq_along(ExonIntron_2)]
    levels(strand(ExonIntron_1)) = levels(strand(ExonIntron_2)) = c("-","+","*") 
    rm(tempCov.antisense)
    cat('Sample- '); cat(i); cat('\n')
    gc()
    
    .resizeCovMat <- function(tempCov, range.index.list, sizeFactor)
    {
      x1 = lapply(range.index.list, function(x) x[[1]])
      x2 = lapply(range.index.list, function(x) x[[2]])
      x3 = lapply(range.index.list, function(x) x[[3]])
      
      out_cov = foreach(idx = list(x1, x2, x3), .combine = c) %dopar%
      {
        cov = sapply(seq_along(tempCov), function(i) { 
          idx_start = idx[[i]][1]; idx_end = idx[[i]][2]
          tempCov[[i]][ idx_start : idx_end ] })
        
        resize_cov = lapply(cov, function(y){ 
          log1p(as.numeric(.covResize(y, resSize)) / sizeFactor)
          })
        resize_cov = Reduce("+", resize_cov) / length(resize_cov) # average over genes
        return(resize_cov)
      }
      .covSmoothen(out_cov, df = 20)
    }
    
    tssCov = .resizeCovMat(tempCov.1, range.index.list, sizeFactor)
    tssCovA = .resizeCovMat(tempCov.1a, range.index.list, sizeFactor) 
    ttsCov = .resizeCovMat(tempCov.2, lapply(range.index.list, function(x) x[4:6]), sizeFactor)
    ttsCovA = .resizeCovMat(tempCov.2a, lapply(range.index.list, function(x) x[4:6]), sizeFactor)
    gc()
    
    cat('Sample: ', i, '\n')
    
    if (sameStrand)
    {
      tssCovArray = c(tssCov, tssCovA)
      ttsCovArray = c(ttsCov, ttsCovA)
    } else {
      tssCovArray = c(tssCovA, tssCov)
      ttsCovArray = c(ttsCovA, ttsCov)
    }
    cov.TEI.list$TSS = cbind(cov.TEI.list$TSS, tssCovArray )
    cov.TEI.list$TTS = cbind(cov.TEI.list$TTS, ttsCovArray )
  }

  return( cov.TEI.list )
}

