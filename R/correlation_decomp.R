`.` <- list # alias to avoid NSE notes in R CMD check




#' @importFrom stats cov.wt
#'
cov_tbl <- function(data, chromosome, signals, rm.boundary = TRUE){
  level <- cov <- weight <- NULL # due to NSE notes in R CMD check
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
  chrcov <- cov.wt(chrmeans[, signals, with = FALSE], wt = chrmeans$weight)$cov[1, 2]

  cov_tbl <- rbind(cov_tbl,
                   data.table(level = chromosome, cov = chrcov))
  return(cov_tbl)
}


# use for testing
# data <- data.table(group = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000))
# data[, y := x + rnorm(4000, mean=0, sd=2)]
# signals <- c("x", "y"); chromosome <- "group"; rm.boundary <- F
# with(data, plot(y ~ x))

#' Correlation decomposition of two signals measured along multiple chromosomes
#'
#' @param data data.table containing columns for variables to be correlated and a column with chromosome id if applicable.
#' data.table should be keyed by chromosome and then position. It is assumed that all positions in the data.table are equally spaced along the chromosome.
#' If not, results will not be valid. To get measurements at equally-spaced locations along chromosomes, you can interpolate your signals using the R function \code{approx()}.
#' @param signals character vector, names of columns to be correlated
#' @param chromosome character string, name of column containing chromosome id. To calculate results per chromosome, use chromosome = NA, and run function separately per chromosome.
#' @param method one of 'pearson' or 'spearman', indicating the type of correlation to calculate from wavelet coefficients. Note 'pearson' wavelet correlations for detail coefficients
#' do not strictly conform to the definition of the Pearson correlation as the product of means is not subtracted, but this gives an unbiased estimate of the wavelet correlation.
#' For scaling coefficients, the usual Pearson formula is used. Note also that the exact decomposition of the total correlation as described in manuscript is only valid using the
#'  'pearson' wavelet and scaling correlations. \cr
#' @param rm.boundary logical, whether to remove boundary coefficients in the calculations. Default is FALSE.
#'
#' @return data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Levels "d?" correspond to wavelet coefficients.
#'  Levels "s?" correspond to scaling coefficients. Note 'pearson' correlations for wavelet coefficients don't subtract off the product of means.
#'  Otherwise 'pearson' correlations for scaling coefficients or at the chromosome-level, and all 'spearman' correlations are calculated in the usual way.
#'  Chromosome-level correlations are weighted by chromosome length.\cr
#'  \tab \cr
#'  \code{cor_n} \tab Correlation computed at the given level \cr
#'  \tab \cr
#'  \code{cor_jack} \tab Bias-corrected correlation from weighted block jackknife procedure. Details of procedure found at # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf \cr
#'  \tab \cr
#'  \code{cor_jack_se} \tab Standard error of correlation computed from a weighted block jackknife across chromosomes. Only returned if number of chromosomes containing the given level is at least 6, otherwise returns NA. Currently only implemented for pearson correlations. \cr
#' }
#'
#' @export
#' @import wCorr
#' @import data.table
#' @importFrom stats cor
#'
#' @examples
#' # Toy example - two white-noise signals "x" and "y" measured
#' # along six 'chromosomes'.
#'
#' library(data.table)
#' data <- data.table(group = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000))
#' data[, y := x + rnorm(4000, mean=0, sd=2)]
#' signals <- c("x", "y"); chromosome <- "group"
#' gnom_cor_decomp(data, signals=signals, chromosome = chromosome)
#'

gnom_cor_decomp <- function(data, signals, chromosome, method = 'pearson', rm.boundary = TRUE){
  level <- weight <- N <- cor_jack <- cor_n <- cor_jack_se <- NULL # due to NSE notes in R CMD check
  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)

  cols <- paste0("coefficient.", signals)

  if (method == 'pearson') {
    # wavelet 'correlations' these are not quite correlations bc we don't subtract off the product of the means
    cor_tbl_detail_prs <- w[grepl("d", level, fixed=T),
                            .(cor_n = mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))), by = level]

    # scaling coefficient correlations do subtract off the product of means
    cor_tbl_smooth_prs <- w[grepl("s", level, fixed=T),
                            .(cor_n = cor(get(cols[1]), get(cols[2]))), by = level]
    cor_tbl <- rbind(cor_tbl_detail_prs, cor_tbl_smooth_prs)

  } else if(method == 'spearman'){
    cor_tbl <- w[, .(cor_n = cor(get(cols[1]), get(cols[2]))), by = level]
  }


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

  if(method == 'pearson'){
    chrcor <- data.table(cor_n = weightedCorr(x = chrmeans[, get(signals[1])],
                                                y = chrmeans[, get(signals[2])],
                                                weights = chrmeans$weight, method = "Pearson"),
                         level = chromosome)

  } else if (method == 'spearman'){
    chrcor <- data.table(cor_n = weightedCorr(x = chrmeans[, get(signals[1])],
                                              y = chrmeans[, get(signals[2])],
                                              weights = chrmeans$weight, method = "Spearman"),
                         level = chromosome)
  }

  cor_tbl <- rbind(cor_tbl, chrcor)

  # ===== weighted jacknife =====
  # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf

  # compute separately for each level, as each level has potentially different number of chromosomes available
  for(j in unique(cor_tbl$level)){

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

    # don't perform jacknife if sample size is at less than 6
    if(g < 6){
      cor_tbl[level == j, c("cor_jack", "cor_jack_se") := .(NA, NA)]

      } else if (g >= 6){
      # leave-1-out estimates will go in this list, first element pearson cors
      # 2nd element spearman cors
      cor_j <- c()

      # loop over chroms and calculate estimate leaving each out in turn
      for(i in 1:length(jchrs)){

        # if detail coefficient
        if(grepl('d', j, fixed=T)){

          cor_j[i] <- w_sub[get(chromosome) != jchrs[i], mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))]
          # smooth coefficient
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
          cor_j[i] <- weightedCorr(x = chrmeans_j[, get(signals[1])],
                                      y = chrmeans_j[, get(signals[2])],
                                      weights = chrmeans_j$weight, method = "Pearson")
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




#' Wavelet partial correlations of two signals measured along multiple chromosomes
#'
#' @param data data.table containing columns for variables to be correlated and a column with chromosome id if applicable.
#' data.table should be keyed by chromosome and then position. It is assumed that all positions in the data.table are equally spaced along the chromosome.
#' If not, results will not be valid. To get measurements at equally-spaced locations along chromosomes, you can interpolate your signals using the R function \code{approx()}.
#' @param signals character vector, names of columns to be correlated. Wavelet coefficients of these variables will be regressed against those of 'z', and then correlations taken.
#' @param z character, name of column containing confounding variable
#' @param chromosome character string, name of column containing chromosome id. To calculate results per chromosome, use chromosome = NA, and run function separately per chromosome.
#' @param method one of 'pearson' or 'spearman', indicating the type of correlation to calculate from wavelet coefficients. Note 'pearson' wavelet correlations for detail coefficients
#' do not strictly conform to the definition of the Pearson correlation as the product of means is not subtracted, but this gives an unbiased estimate of the wavelet correlation.
#' For scaling coefficients, the usual Pearson formula is used. Note also that the exact decomposition of the total correlation as described in manuscript is only valid using the
#'  'pearson' wavelet and scaling correlations. \cr
#' @param rm.boundary logical, whether to remove boundary coefficients in the calculations. Default is FALSE.
#'
#' @return data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Levels "d?" correspond to wavelet coefficients.
#'  Levels "s?" correspond to scaling coefficients. Note 'pearson' correlations for wavelet coefficients don't subtract off the product of means.
#'  Otherwise 'pearson' correlations for scaling coefficients or at the chromosome-level, and all 'spearman' correlations are calculated in the usual way.
#'  Chromosome-level correlations are weighted by chromosome length.\cr
#'  \tab \cr
#'  \code{cor_n} \tab Correlation computed at the given level \cr
#'  \tab \cr
#'  \code{cor_jack} \tab Bias-corrected correlation from weighted block jackknife procedure. Details of procedure found at # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf \cr
#'  \tab \cr
#'  \code{cor_jack_se} \tab Standard error of correlation computed from a weighted block jackknife across chromosomes. Only returned if number of chromosomes containing the given level is at least 6, otherwise returns NA. \cr
#' }
#'
#' @export
#' @import wCorr
#' @import data.table
#' @importFrom stats cor
#'
#' @examples
#' # Toy example - two white-noise signals "x" and "y" measured
#' # along six 'chromosomes'.
#'
#' library(data.table)
#' data <- data.table(group = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000))
#' data[, y := x + rnorm(4000, mean=0, sd=2)]
#' data[, z := y + rnorm(4000, mean=0, sd=2)]
#' signals <- c("x", "y"); chromosome <- "group"
#' z <- "z"
#' gnom_partial_cor_decomp(data, signals=signals, z = z, chromosome = chromosome)
#'
gnom_partial_cor_decomp <- function(data, signals, z, chromosome, method = 'pearson', rm.boundary = TRUE){
  level <- weight <- N <- cor_jack <- cor_n <- cor_jack_se <- NULL # due to NSE notes in R CMD check
  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = c(signals, z), rm.boundary = rm.boundary)

  cols1 <- paste0("coefficient.", c(signals, z))

  w[, residuals.x := residuals(lm(get(cols1[1]) ~ get(cols1[3]))), by = level]
  w[, residuals.y := residuals(lm(get(cols1[2]) ~ get(cols1[3]))), by = level]

  cols <- paste0("residuals.", c("x", "y", "z"))

  if (method == 'pearson') {
    # wavelet 'correlations' these are not quite correlations bc we don't subtract off the product of the means

    cor_tbl_detail_prs <- w[grepl("d", level, fixed=T),
                            .(cor_n = mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))), by = level]

    # scaling coefficient correlations do subtract off the product of means
    cor_tbl_smooth_prs <- w[grepl("s", level, fixed=T),
                            .(cor_n = cor(get(cols[1]), get(cols[2]))), by = level]
    cor_tbl <- rbind(cor_tbl_detail_prs, cor_tbl_smooth_prs)

  } else if(method == 'spearman'){
    cor_tbl <- w[, .(cor_n = cor(get(cols[1]), get(cols[2]))), by = level]
  }


  if(is.na(chromosome)){
    return(cor_tbl)
  }

  # chromosome-scale correlation:
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.", signals))
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = c(signals,z)]

  # chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)

  chrmeans[, residuals.x := residuals(lm(get(signals[1])~get(z), weights = weight))]
  chrmeans[, residuals.y := residuals(lm(get(signals[2])~get(z), weights = weight))]

  # 4. weighted chromosome-scale correlation

  if(method == 'pearson'){
    # since residuals already weighted, just compute regular correlation
    chrcor <- chrmeans[, .(cor_n = cor(residuals.x, residuals.y), level = chromosome)]

  } else if (method == 'spearman'){
    chrcor <- chrmeans[, .(cor_n = cor(residuals.x, residuals.y, method = 'spearman'), level = chromosome)]

  }

  cor_tbl <- rbind(cor_tbl, chrcor)

  # ===== weighted jacknife =====
  # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf

  # compute separately for each level, as each level has potentially different number of chromosomes available
  for(j in unique(cor_tbl$level)){

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

    # don't perform jacknife if sample size is at less than 6
    if(g < 6){
      cor_tbl[level == j, c("cor_jack", "cor_jack_se") := .(NA, NA)]

    } else if (g >= 6){
      # leave-1-out estimates will go in this list, first element pearson cors
      # 2nd element spearman cors
      cor_j <- c()

      # loop over chroms and calculate estimate leaving each out in turn
      for(i in 1:length(jchrs)){

        # if detail coefficient
        if(grepl('d', j, fixed=T)){

          cor_j[i] <- w_sub[get(chromosome) != jchrs[i], mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))]
          # smooth coefficient
        } else if(grepl('s', j, fixed=T)){

          cor_j[i] <- w_sub[get(chromosome) != jchrs[i], cor(get(cols[1]), get(cols[2]))]

        } else if(j == chromosome){

          # chromosome-scale correlation:
          totalmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), .SDcols = signals]
          setnames(totalmeans_j, signals, paste0("totalmean.", signals))
          chrmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), by = chromosome, .SDcols = c(signals, z)]

          # chromosome lengths as weights
          chrlen_j <- data[get(chromosome) != jchrs[i], .(weight = .N), by = chromosome]
          chrlen_j[, weight := weight/sum(weight)]
          chrmeans_j <- merge(chrmeans_j, chrlen_j)

          #  weighted chromosome-scale correlation
          chrmeans_j[, residuals.x := residuals(lm(get(signals[1])~get(z), weights = weight))]
          chrmeans_j[, residuals.y := residuals(lm(get(signals[2])~get(z), weights = weight))]

          cor_j[i] <- chrmeans[, cor(residuals.x, residuals.y)]

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
#' @param rm.boundary logical, whether to remove boundary coefficients (default TRUE)
#'
#' @return if avg.over.chroms = TRUE, returns a data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Values "dX"
#' correspond to wavelet coefficients, and corresponding variances are associated
#'  with changes in localized average signal values over scale 2^(X-1). Values "sX" correspond to
#'  scaling coefficients, and corresponding variances are associated with localized averages
#'  on a scale of 2^X.  \cr
#'  \tab \cr
#'  \code{rsqrd_n} \tab R squared from a linear model using the specified formula.
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
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @export
#'
#' @examples
#' # sample data with 20 chroms, where z is correlated with both x and y separately
#'
#' library(data.table)
#'
#' chrlens <- sample(x = c(100, 200, 300, 400, 500), size = 20, replace=TRUE)
#' chr_id <- NULL
#' for(i in 1:20){
#'   chr_id <- c(chr_id, rep(i, chrlens[i]))
#' }
#' data <- data.table(chrom = chr_id, x = rnorm(length(chr_id)), y = rnorm(length(chr_id)))
#' data[, z := 4*x - 2*y + rnorm(length(chr_id), sd = 5)]
#'
#' modwt_lm_rsqrd(data=data, yvar = "z", xvars = c("x", "y"), chromosome = "chrom")
#'
modwt_lm_rsqrd <- function(data, yvar, xvars, chromosome, rm.boundary = TRUE){
  level <- N <- rsqrd_n <- rsqrd_jack <- rsqrd_jack_se <- NULL # due to NSE notes in R CMD check
  # ===== lm of wavelet coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = c(xvars,yvar), rm.boundary = rm.boundary)

  # formula for lm
  f_detail <- as.formula(paste(paste0('coefficient.',yvar), paste(paste0('coefficient.',xvars), collapse=" + "), sep=" ~ "))
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

    if(g < 6){
      rsqrd[level == j, c("rsqrd_jack", "rsqrd_jack_se") := .(NA, NA)]
    } else if (g >=6){

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
