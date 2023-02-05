
cov_tbl <- function(data, chromosome, signals, rm.boundary = FALSE){
  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  cols <- paste0("coefficient.", signals)

  # wavelet covariances
  cov_tbl_detail <- w[grepl("d", level, fixed=T), .(cov = mean(get(cols[1])*get(cols[2]))), by = level]
  cov_tbl_smooth <- w[grepl("s", level, fixed=T), .(cov = cov(get(cols[1]), get(cols[2]))), by = level]
  cov_tbl <- rbind(cov_tbl_detail, cov_tbl_smooth)

  if(is.na(chromosome)){
    return(cov_tbl)
  }

  # chromosome-scale covariance:
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.", signals))
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)

  # 4. weighted chromosome-scale covariance
  chrcov <- cov.wt(chrmeans[, ..signals ], wt = chrmeans$weight)$cov[1, 2]

  cov_tbl <- rbind(cov_tbl,
                   data.table(level = chromosome, cov = chrcov))
  return(cov_tbl)
}


# use for testing
# data <- data.table(group = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000))
# data[, y := x + rnorm(4000, mean=0, sd=2)]
# signals <- c("x", "y")
# with(data, plot(y ~ x))

gnom_cor_decomp <- function(data, chromosome, signals, rm.boundary = FALSE){

  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)

  cols <- paste0("coefficient.", signals)

  # wavelet 'correlations' these are not quite correlations bc we don't subtract off the product of the means
  cor_tbl_detail <- w[grepl("d", level, fixed=T), .(cor_n = mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))), by = level]

  # scaling coefficient correlations do subtract off the product of means
  cor_tbl_smooth <- w[grepl("s", level, fixed=T), .(cor_n = cor(get(cols[1]), get(cols[2]))), by = level]
  cor_tbl <- rbind(cor_tbl_detail, cor_tbl_smooth)

  if(is.na(chromosome)){
    return(cor_tbl)
  }

  # chromosome-scale correlation:
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.", signals))
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)

  # 4. weighted chromosome-scale correlation
  chrcor <- cov.wt(chrmeans[, ..signals ], wt = chrmeans$weight, cor = TRUE)$cor[1, 2]

  cor_tbl <- rbind(cor_tbl,
                  data.table(level = chromosome, cor_n = chrcor))

  # ===== weighted jacknife =====
  # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf

  # compute separately for each level, as each level has potentially different number of chromosomes
  for(j in cor_tbl$level){

    if(j != chromosome){
      # which chromosomes have coefficients at this level
      w_sub <- w[level == j]
      n <- nrow(w_sub)
      jchrs <- w_sub[, unique(get(chromosome))]

      # number of observations for each chromosome where this level is present
      j_weights <- w_sub[, .N, by = chromosome][, N]

      # subset data to only these chromosomes
      d_sub <- data[get(chromosome) %in% jchrs]

      # how many chromosomes in our sample for this level?
      g <- length(jchrs)

    } else if (j == chromosome){
      n <- nrow(data)
      jchrs <- data[, unique(get(chromosome))]
      j_weights <- data[, .N, by = chromosome][, N]
      g <- length(jchrs)
    }


    if(g < 4){
      cor_tbl[level == j, c("cor_jack", "cor_jack_se") := .(NA, NA)]
    } else if (g >=4){
      # leave-1-out estimates will go in this vector
      cor_j <- vector()

      # loop over chroms and calculate estimate leaving each out in turn
      for(i in 1:length(jchrs)){
        # if detail coefficient
        if(grepl('d', j, fixed=T)){
          cor_j[i] <- w_sub[get(chromosome) != jchrs[i], mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))]
        } else if(grepl('s', j, fixed=T)){
          cor_j[i] <- w_sub[get(chromosome) != jchrs[i], cor(get(cols[1]), get(cols[2]))]
        } else if(j == chromosome){

          # chromosome-scale correlation:
          totalmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), .SDcols = signals]
          setnames(totalmeans_j, signals, paste0("totalmean.", signals))
          chrmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), by = chromosome, .SDcols = signals]

          # chromosome lengths as weights
          chrlen_j <- data[get(chromosome) != jchrs[i], .(weight = .N), by = chromosome]
          chrlen_j[, weight := weight/sum(weight)]
          chrmeans_j <- merge(chrmeans_j, chrlen_j)

          #  weighted chromosome-scale correlation
          cor_j[i] <- cov.wt(chrmeans_j[, ..signals ], wt = chrmeans_j$weight, cor = TRUE)$cor[1, 2]

        }
      }

      # compute jacknife weighted estimate of the correlation
      cor_tbl[level==j, cor_jack := g*cor_n - sum(((n-j_weights)*cor_j)/n)]

      # compute jacknife standard error (formula in link gives jacknife variance of the estimator, thus we just take square root for se)
      h_j <- n/j_weights
      tau_j <- h_j*cor_tbl[level==j, cor_n] - (h_j-1)*cor_j

      cor_tbl[level==j, cor_jack_se := sqrt((1/g)*sum((tau_j - cor_j)^2/(h_j-1)))][]
    }
  }

  return(cor_tbl)
}


# make sample data set with 20 chroms, where z is correlated with both x and y separately
# chrlens <- sample(x = c(100, 200, 300, 400, 500), size = 20, replace=T)
# chr_id <- NULL
# for(i in 1:20){
#   chr_id <- c(chr_id, rep(i, chrlens[i]))
# }
# data <- data.table(chrom = chr_id, x = rnorm(length(chr_id)), y = rnorm(length(chr_id)))
# data[, z := 4*x - 2*y + rnorm(length(chr_id), sd = 5)]
# chromosome <- 'chrom'
# rm.boundary <- F
# xvars <- c("x", "y")
# yvar <- "z"


#' R squared from linear model of wavelet coefficients
#'
#' Estimates r squared from a linear model of wavelet coefficients.
#' If separate chromosomes are present, also estimates r squared from linear model of chromosome means,
#' weighting by chromosome length. Weighted block jacknife with chromosomes as blocks
#' is performed to return bias-corrected estimates and standard errors of the bias-corrected estimate.
#'
#' @param data data.frame containing columns for chromosome id and signal measurements.
#' It is assumed that the rows contain datapoints that are evenly spaced along the sequence.
#' If not, you will need to interpolate your signal to a grid of evenly spaced points.
#' @param yvar character string, the response variable in the linear model
#' @param xvars character vector, the predictor variables in the linear model. Assumed no interaction
#' @param chromosome character string, the name of the column containing chromosome id
#' @param rm.boundary logical, whether to remove boundary coefficients (default FALSE)
#'
#' @return if avg.over.chroms = TRUE, returns a data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Values "dX"
#' correspond to wavelet coefficients, and corresponding variances are associated
#'  with changes in localized average signal values over scale 2^(X-1). Values "sX" correspond to
#'  scaling coefficients, and corresponding variances are associated with localized averages
#'  on a scale of 2^X.  \cr
#'  \tab \cr
#'  \code{rsqrd_n} \tab R squared from a linear model using the specified formula.
#'  For wavelet coefficients, the regression is forced through the origin.
#'  For scaling coefficients and chromosome means, an intercept is fit.
#'  For the chromosome level, the regression of means is weighted by chromosome length.
#'  \cr
#'  \tab \cr
#'  \code{rsqrd_jack} \tab The bias-corrected estimate of r squared from a weighted block jacknife with chromosomes as blocks.
#'  Not reported if the given level contains < 4 blocks.
#'   \cr
#'   \tab \cr
#'  \code{rsqrd_jack_se} \tab Standard error of the jacknife estimate.
#' }
#'
#'
#' @export
#'
#' @examples
#' # sample data with 20 chroms, where z is correlated with both x and y separately
#'
#' chrlens <- sample(x = c(100, 200, 300, 400, 500), size = 20, replace=T)
#' chr_id <- NULL
#' for(i in 1:20){
#'   chr_id <- c(chr_id, rep(i, chrlens[i]))
#' }
#' data <- data.table(chrom = chr_id, x = rnorm(length(chr_id)), y = rnorm(length(chr_id)))
#' data[, z := 4*x - 2*y + rnorm(length(chr_id), sd = 5)]
#'
#' modwt_lm_rsqrd(data=data, yvar = "z", xvars = c("x", "y"), chromosome = "chrom")
#'
modwt_lm_rsqrd <- function(data, yvar, xvars, chromosome, rm.boundary = FALSE){

  # ===== lm of wavelet coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = c(xvars,yvar), rm.boundary = rm.boundary)

  # formula for lm, intercept forced through zero
  f_detail <- as.formula(paste(paste0('coefficient.',yvar), paste(c('0',paste0('coefficient.',xvars)), collapse=" + "), sep=" ~ "))
  f_smooth <- as.formula(paste(paste0('coefficient.',yvar), paste(paste0('coefficient.',xvars), collapse=" + "), sep=" ~ "))
  f_chr <- as.formula(paste(yvar, paste(xvars, collapse=" + "), sep =" ~ "))

  rsqrd_d <- w[grepl('d',level,fixed=T), .(rsqrd_n = summary(lm(f_detail, data = .SD))$r.squared), by = level]
  rsqrd_s <- w[grepl('s',level,fixed=T), .(rsqrd_n = summary(lm(f_smooth, data = .SD))$r.squared), by = level]
  rsqrd <- rbind(rsqrd_d, rsqrd_s)

  if(is.na(chromosome)){
    return(rsqrd[])
  }

  # rsquared from weighted lm of chromosome means
  chrmeans <- data[, lapply(.SD, mean), .SDcols = c(xvars, yvar), by = chromosome]
  chr_weights <- data[, .(weight = .N), by = chromosome]
  chrmeans <- merge(chrmeans, chr_weights)
  rsqrd_chr <- summary(lm(f_chr, weights = chrmeans$weight, data = chrmeans))$r.squared
  rsqrd <- rbind(rsqrd, data.table(level = chromosome, rsqrd_n = rsqrd_chr))

  # weighted jacknife for wavelet coefficients
  for(j in rsqrd$level){
    if(j != chromosome){
      # which chromosomes have coefficients at this level
      w_sub <- w[level == j]
      n <- nrow(w_sub)
      jchrs <- w_sub[, unique(get(chromosome))]

      # number of observations for each chromosome where this level is present
      j_weights <- w_sub[, .N, by = chromosome][, N]

      # subset data to only these chromosomes
      d_sub <- data[get(chromosome) %in% jchrs]

      # how many chromosomes in our sample for this level?
      g <- length(jchrs)

    } else if (j == chromosome){
      n <- nrow(data)
      jchrs <- data[, unique(get(chromosome))]
      j_weights <- data[, .N, by = chromosome][, N]
      g <- length(jchrs)
    }

    if(g < 4){
      rsqrd[level == j, c("rsqrd_jack", "rsqrd_jack_se") := .(NA, NA)]
    } else if (g >=4){

      # leave-1-out estimates will go in this vector
      rsqrd_j <- vector()

      # loop over chroms and calculate estimate leaving each out in turn
      for(i in 1:length(jchrs)){

        # if detail coefficient
        if(j != chromosome & grepl('d', j, fixed=T)){
          rsqrd_j[i] <- w_sub[get(chromosome) != jchrs[i], summary(lm(f_detail, data = .SD))$r.squared]
        } else if(j != chromosome & grepl('s', j, fixed=T)){
          rsqrd_j[i] <- w_sub[get(chromosome) != jchrs[i], summary(lm(f_smooth, data = .SD))$r.squared]
        } else if(j == chromosome){

          # rsquared from weighted lm of chromosome means
          chrmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), .SDcols = c(xvars, yvar), by = chromosome]
          chr_weights_j <- data[get(chromosome) != jchrs[i], .(weight = .N), by = chromosome]
          chrmeans_j <- merge(chrmeans_j, chr_weights_j)
          rsqrd_j[i] <- summary(lm(f_chr, weights = chrmeans_j$weight, data = chrmeans_j))$r.squared
        }
      }

      # compute jacknife weighted estimate of the correlation
      rsqrd[level==j, rsqrd_jack := g*rsqrd_n - sum(((n-j_weights)*rsqrd_j)/n)][]

      # compute jacknife standard error (formula in link gives jacknife variance of the estimator, thus we just take square root for se)
      h_j <- n/j_weights
      tau_j <- h_j*rsqrd[level==j, rsqrd_n] - (h_j-1)*rsqrd_j
      rsqrd[level==j, rsqrd_jack_se := sqrt((1/g)*sum((tau_j - rsqrd_j)^2/(h_j-1)))][]
    }
  }
  return(rsqrd[])
}
