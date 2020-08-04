## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## predict copy number
weightTransform <- function( W ) {
  out_W <- 2^W - 0.001
  out_W[out_W < 0] <- 0
  out_W
}

# Constants for RNA copy number transformation
# 1) correct spike-in weight due to the different mole number design to adapt to 
# the final mixed concentration 
#     before      now
# Sp2      1        1
# Sp4      1      0.1
# Sp5      1        1
# Sp8      1      0.1
# Sp9      1     0.01
# Sp12     1     0.01
# sum      6     2.22
norm_factor <- 6 / 2.22 

# 2) Transform 0.4 ng spike-in RNA mixture per million cell to molecular number per cell:
# ( 0.4 * 1e-9 (g) / 322 (nt molecular mass) * 6.02e+23 ) / 1000 (to RPK average length) / 1e+6 (million cells) / norm_factor 
transform_facter <- 0.4e-9 / 322 / 1000 * 6.02e+23 * 1e-6 / norm_factor

## predict labeled rate and half-life
PredictRate <- function(beta, LRNA, FRNA) {
  rate_term <- cbind(1, FRNA, LRNA) %*% beta
  rate <- 2^(rate_term) - 0.001
  rate[is.na(rate) | rate < 0] <- 0
  rate
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# geomatric mean by row
geoMeans <- function(x) exp( log(matrixStats::rowProds(x)) / ncol(x))

sampleMatMeans <- function( sample_counts, FUN ) {
  mat <- matrix(sample_counts, ncol = 9)
  cbind(sample_counts[, 1:2] %>% FUN,
        sample_counts[, 3],
        sample_counts[, 4:5] %>% FUN,
        sample_counts[, 6:7] %>% FUN,
        sample_counts[, 8:9] %>% FUN
  )
}
