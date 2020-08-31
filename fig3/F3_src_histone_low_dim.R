source("../util/getCoverage.R")

# functions ----------------------------------------------------------------------------------------

trim_quantile <- function(x, q = 0.999)
{
  x[is.na(x)] <- 0
  x[x > quantile(x, q)] <- quantile(x, q)
  x[x < quantile(x, 1 - q)] <- quantile(x, 1 - q)
  return(x)
}


joint_Projection <- function(con1_mat, con2_mat)
{
  # This function projects two matrices on the average matrix eigen vector,
  # to compare difference with the first two principle components
  # no scaling applied since features are inherently different distributed ChIP signals
  mean_mat <- (con1_mat + con2_mat) / 2
  
  mean_eigen_vector <- mean_mat %>% cov() %>% eigen() %>% "$"("vectors")

  con1_pc <- as.matrix(sweep(con1_mat, 2, colMeans(con1_mat), "-")) %*% mean_eigen_vector
  con2_pc <- as.matrix(sweep(con2_mat, 2, colMeans(con2_mat), "-")) %*% mean_eigen_vector

  list(PCs = data.frame(PC1.1 = con1_pc[, 1],
                        PC2.1 = con1_pc[, 2],
                        PC1.2 = con2_pc[, 1],
                        PC2.2 = con2_pc[, 2]),
       mean_eigen_vector = mean_eigen_vector)
}

grid_arrow <- function(dat, grid_n) {
  # from dat[, 1:2] to dat[, 3:4]
  grid_step_x = diff(range(dat[,1])) / (grid_n + 2)
  grid_step_y = diff(range(dat[,2])) / (grid_n + 2)
  grid_dat <- expand.grid(seq_len(grid_n), seq_len(grid_n))
  grid_dat$x1 <- seq(min(dat[,1]), max(dat[,1]), diff(range(dat[,1])) / (grid_n + 2))[grid_dat$Var1 + 1]
  grid_dat$y1 <- seq(min(dat[,2]), max(dat[,2]), diff(range(dat[,2])) / (grid_n + 2))[grid_dat$Var2 + 1]
  grid_dat$y2 = grid_dat$x2 = NA
  for (i in seq_len(nrow(grid_dat))) {
    idx <- dat[,3] > (min(dat[,3]) + grid_step_x * (grid_dat[i, 1]-1)) &
      (dat[,3] < (min(dat[,3]) + grid_step_x * (grid_dat[i, 1]+1))) &
      (dat[,4] > (min(dat[,4]) + grid_step_y * (grid_dat[i, 2]-1))) &
      (dat[,4] < (min(dat[,4]) + grid_step_y * (grid_dat[i, 2]+1)))
    grid_dat$x2[i] <- mean(dat[,3][idx])
    grid_dat$y2[i] <- mean(dat[,4][idx])
  }
  grid_dat
}

plot_PCA_shifts <- function(mat_1, mat_2,
                            directions = c(1, 1),
                            gene_idx = NULL,
                            idx_color = "black",
                            grid_n = 20,
                            eigen_vector_color,
                            point_scale = NULL,
                            point_scale_color,
                            scale_name = "log2FC")
{ # Args:
  # mat_1, mat_2: matrices of observations to be compared
  # gene_idx: gene indexes to be pinpointed
  # eigen_vector_color: temporary colors with names matching with Marks
  # grid_n: number of arrows on each PC dimension
  # point_scale: log2FC for PCs points color gradient
  # scale_name: legend tile

  # first 2 PCs
  pc_dat <- joint_Projection(mat_1, mat_2) %>% "$"(PCs) %>%
    as.matrix() %*% diag(rep(directions, 2)) %>%
    as.data.frame()
  colnames(pc_dat) <- c("PC1", "PC2", "PC1_v", "PC2_v")
  
  # feature eigen vectors
  eigen_vector <- joint_Projection(mat_1, mat_2)$mean_eigen_vector * 2
  ev2_dat <- data.frame(V1 = eigen_vector[, 1] * directions[1],
                        V2 = eigen_vector[, 2] * directions[2],
                        V1_0 = 0,
                        V2_0 = 0,
                        Marks = colnames(mat_1))
  ev2_dat$Marks = gsub("me", "m", ev2_dat$Marks)

  # append plot elements
  p <- ggplot() 
  
  if (!is.null(point_scale)) {
    pc_dat <- cbind(pc_dat, point_scale = point_scale)
    p <- p + geom_point(data = pc_dat[!is.na(pc_dat$point_scale), ],
                        aes(x = PC1, y = PC2,
                            color = cut(point_scale, breaks = c(-Inf, -1, 1, Inf)) ),
                        size = 1) +
      scale_color_manual(values = point_scale_color, name = scale_name)
  }
  
  if (grid_n > 0) {
    # grid of shifting
    grid_dat <- grid_arrow(pc_dat, grid_n)
    
    if (!exists("eigen_vector_color"))
      eigen_vector_color <- setNames(RColorBrewer::brewer.pal(9, "Set1"),
                                     c("H3K27ac", "H3K27m2", "H3K27m3", "H3K4m3", "H4K12ac", "H4K5ac",
                                       "Nanog", "Pou5F1", "Sox2"))
    
    p <- p + geom_segment(data = grid_dat,
                          mapping = aes(x = x1, y = y1, xend = x2, yend = y2),
                          arrow = arrow(type="open", angle=20, length = unit(x = 0.13, 'cm')),
                          size = 0.25, color = add.alpha("darkblue", 0.5)) +
      geom_hline(yintercept=0, linetype="dashed",
                 color = "grey", size= .5) +
      geom_vline(xintercept=0, linetype="dashed",
                 color = "grey", size= .5) +
      geom_segment(data = ev2_dat,
                   mapping = aes(x = V1_0, y = V2_0, xend = V1, yend = V2),
                   arrow = arrow(type="open", angle=30, length = unit(x = 0.3, 'cm')),
                   col = eigen_vector_color[ev2_dat$Marks],
                   size = 1) +
      scale_color_manual(values = eigen_vector_color[as.character(ev2_dat$Marks)])
  } 
  if (!is.null(gene_idx)) {
    p <- p + geom_point(data = pc_dat[gene_idx,], aes(x = PC1, y = PC2),
                        size = 1.5, color = add.alpha("red", 0.5), pch = 17)
  }
  p + theme_setting + 
      theme(panel.border = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size = 0.5)))
}




plot_PCA_contrast <- function(mat_1, mat_2, gene_idx) {

  pc = prcomp(mat_2 - mat_1)

  loadings <- data.frame(pc$rotation[, 1:2],
                         Marks = colnames(mat_1))

  # set color only for the current plot
  loadings$Marks = gsub("me", "m", loadings$Marks)
  eigen_vector_color <- setNames(RColorBrewer::brewer.pal(9, "Set1"),
                        c("H3K27ac", "H3K27m2", "H3K27m3", "H3K4m3", "H4K12ac", "H4K5ac",
                          "Nanog", "Pou5F1", "Sox2"))

  ggplot(data.frame(pc$x), aes(x = PC1, y = PC2)) +
    geom_hex() +
    geom_hline(yintercept=0, linetype="dashed",
               color = "grey", size= .5) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey", size= .5) +
    geom_point(data = data.frame(pc$x)[gene_idx,],
               aes(x = PC1, y = PC2),
               size = 0.4, color = rgb(0,0,0, 1)) +
    scale_fill_viridis(option = "E", direction = -1) +
    labs(fill = "Genes") +
    geom_segment(data = loadings,
                 mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Marks),
                 arrow = arrow(type="open", angle=20, length = unit(x = 0.2, 'cm')),
                 size = 1) +
    scale_color_manual(values = eigen_vector_color[as.character(loadings$Marks)])
}


# load MINUTE ChIP samples ------------------------------------------------------------------------------------------------
tss.gr <- importRanges("../data/tss.attributes.gff3")
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(tss.gr$gene_name),
                       keytype = "GENENAME",
                       columns = c("ENTREZID", "GENEID"))
tss.gr$gene_id <- res$GENEID[match(tss.gr$gene_name, res$GENENAME)]
tss.gr <- tss.gr[!is.na(tss.gr$gene_id)]

bw_files <- list.files('/mnt/E0767589767560E8/UPPMAX/SL2i',
                       pattern = 'H3|H4.*bw', full.names = T)

design_mat <- data.frame(condition = gsub("(.*_)(.*)_(CTR|INH).*", "\\2", bw_files),
                         mark = gsub(".*SL2i/(H.*)_(SL|2i).*", "\\1", bw_files),
                         treat = gsub("(.*_)(.*)_.*", "\\2", bw_files))

# get histone mark counts
if (file.exists("data/tss_histone_counts.RData")) {
  # load data
  tss_histone_counts <- readRDS("data/tss_histone_counts.RData")
  design_mat <- readRDS("data/MINUTE_design_mat.RData")
  tss_histone_counts2 <- readRDS("data/tss_histone_counts2.RData")
  design_mat2 <- readRDS("data/Galonska_design_mat2.RData")
  tss_histone_counts3 <- readRDS("data/tss_histone_counts3.RData")
  design_mat3 <- readRDS("data/Marks_design_mat3.RData")

  # Select marks to compare and scale counts

  # SL vs 2i: H3K27m3, H3K4m3, H3K27m2
  mat_SL <- tss_histone_counts[, design_mat$treat == "CTR" &
                                 design_mat$condition == "SL" &
                                 design_mat$mark %in% c("H3K27m2", "H3K27m3", "H3K4m3")] %>%
    as.data.frame() %>% 
    # mutate_all(funs(.-mean(.))) %>% # substracting colMeans doesn't affect PCA
    `colnames<-`(design_mat$mark[design_mat$treat == "CTR" &
                                   design_mat$condition == "SL" &
                                   design_mat$mark %in% c("H3K27m2", "H3K27m3", "H3K4m3")])

  mat_2i <- tss_histone_counts[, design_mat$treat == "CTR" &
                                 design_mat$condition == "2i" ] %>%
    as.data.frame() %>%
    # mutate_all(funs(.-mean(.))) %>%
    `colnames<-`(design_mat$mark[design_mat$treat == "CTR" &
                                   design_mat$condition == "2i"])

  # SL vs 2i: H3K27m3, H3K4m3, Nanog, Pou5F1, Sox2
  mat_SL.2 <- tss_histone_counts2[, design_mat2$condition == "Serum"] %>%
    as.data.frame() %>% 
    # mutate_all(funs(.-mean(.))) %>%
    `colnames<-`(design_mat2$mark[design_mat2$condition == "Serum"])

  mat_2i.2 <- tss_histone_counts2[, design_mat2$condition == "2i"] %>%
    as.data.frame() %>% 
    # mutate_all(funs(.-mean(.))) %>%
    `colnames<-`(design_mat2$mark[design_mat2$condition == "2i"])

  # SL vs mTORi: H3K27ac, H3K27m2, H3K27m3, H3K4m3, H4K12ac, H4K5ac
  mat_SL_CTR <- tss_histone_counts[, design_mat$treat == "CTR" &
                                     design_mat$condition == "SL"] %>%
    as.data.frame() %>% 
    # mutate_all(funs(.-mean(.))) %>%
    `colnames<-`(design_mat$mark[design_mat$treat == "CTR" &
                                   design_mat$condition == "SL"])

  mat_SL_INH <- tss_histone_counts[, design_mat$treat == "INH" &
                                     design_mat$condition == "SL"] %>%
    as.data.frame() %>% 
    # mutate_all(funs(.-mean(.))) %>%
    `colnames<-`(design_mat$mark[design_mat$treat == "INH" &
                                   design_mat$condition == "SL"])

  } else {
  tss_histone_counts <- .countBW(bw_files = bw_files,
                                 intervals = tss.gr + 1000,
                                 fast = T) %>%
    "*"(1000 / 50) %>% # convert to reads number
    log1p() %>%
    apply(2, trim_quantile, q = 0.999) %>%
    as.data.frame()

  saveRDS(tss_histone_counts, "data/tss_histone_counts.RData")
  saveRDS(design_mat, "data/MINUTE_design_mat.RData")

  # load 2015_Galonska samples--------------------------------------------------------------------------------------------------
  bw_files2 <- list.files('/mnt/0E471D453D8EE463/GEO_bw/2015_Galonska',
                          pattern = '(Serum|3p_2i|Pou5F1_2i).*bw', full.names = T)

  # get histone mark counts
  tss_histone_counts2 <- .countBW(bw_files = bw_files2,
                                  intervals = tss.gr + 1000,
                                  fast = T) %>%
    "*"(1000 / 50) %>%
    log1p() %>%
    apply(2, trim_quantile, q = 0.999) %>%
    as.data.frame()

  tss_histone_counts2$'2015_Galonska_mES_H3K4me3_2i.mm9.bw' <-
    tss_histone_counts2[, 3:4] %>% exp %>% rowMeans %>% log
  tss_histone_counts2$'2015_Galonska_mES_H3K4me3_Serum.mm9.bw' <-
    tss_histone_counts2[, 5:6] %>% exp %>% rowMeans %>% log
  tss_histone_counts2$'2015_Galonska_mES_Nanog_2i.mm9.bw' <-
    tss_histone_counts2[, 7:8] %>% exp %>% rowMeans %>% log
  tss_histone_counts2$'2015_Galonska_mES_Nanog_Serum.mm9.bw' <-
    tss_histone_counts2[, 9:10] %>% exp %>% rowMeans %>% log
  tss_histone_counts2$'2015_Galonska_mES_Sox2_Serum.mm9.bw' <-
    tss_histone_counts2[, 14:15] %>% exp %>% rowMeans %>% log
  tss_histone_counts2$'2015_Galonska_mES_Sox2_2i.mm9.bw' <-
    tss_histone_counts2[, 13]
  tss_histone_counts2$'2015_Galonska_mES_H3K27me3_2i.mm9.bw' <-
    tss_histone_counts2$'2015_Galonska_mES_H3K27me3_3p_2i.mm9.bw'

  tss_histone_counts2 <- tss_histone_counts2[, !grepl("r1|r2|3p", colnames(tss_histone_counts2))]
  tss_histone_counts2 <- tss_histone_counts2[, order(colnames(tss_histone_counts2))]

  design_mat2 <- data.frame(condition = gsub("(.*_)(.*).mm.*",
                                             "\\2", colnames(tss_histone_counts2)),
                            mark = gsub("2015_Galonska_mES_(.*)_(3p|Serum|2i).*",
                                        "\\1", colnames(tss_histone_counts2)))

  saveRDS(tss_histone_counts2, "data/tss_histone_counts2.RData")
  saveRDS(design_mat2, "data/Galonska_design_mat2.RData")

  # load 2012_Marks samples -----------------------------------------------------------------------------------------------------
  bw_files3 <- list.files('/mnt/0E471D453D8EE463/GEO_bw/2012_Marks',
                          pattern = 'E14.*bw', full.names = T)

  # get histone mark counts
  tss_histone_counts3 <- .countBW(bw_files = bw_files3,
                                  intervals = tss.gr + 1000,
                                  fast = T) %>%
    "*"(1000 / 50) %>%
    log1p() %>%
    apply(2, trim_quantile, q = 0.999) %>%
    as.data.frame()

  design_mat3 <- data.frame(condition = gsub("(.*-)(.*)_.*",
                                             "\\2", colnames(tss_histone_counts3)),
                            mark = gsub(".*_(.*).mm.*",
                                        "\\1", colnames(tss_histone_counts3)))

  saveRDS(tss_histone_counts3, "data/tss_histone_counts3.RData")
  saveRDS(design_mat3, "data/Marks_design_mat3.RData")
}

#-----------------------------------------------------------------------------------------------------
