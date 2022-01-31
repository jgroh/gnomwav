sample_var <- function(u){
  return(mean(u^2, na.rm=T) - mean(u, na.rm=T)^2)
}

#' Genome-wide variance decomposition
#'
#' Performs a scale-based decomposition of the variance in signal measurements along a genome
#' containing multiple chromosomes. Wavelet variances are obtained using the Maximal Overlap
#' Discrete Wavelet Transform with Haar wavelets. By default boundary coefficients are
#' removed, which provides an unbiased estimate of the wavelet variance. A weighted average
#' of wavelet variances is taken over chromosmes, and chromosome-level variance is also
#' calculated.
#'
#' @param x data.frame containing columns for chromosome id and signal measurements
#' @param chromosome character string, name of column containing chromosome id
#' @param signals character vector, names of columns containing the signals for which
#' variance decomposition is desired
#' @param rm.boundary logical, whether to remove boundary coefficients
#' @param avg.over.chroms logical, if FALSE, returns variance decomposition per chromosome
#' if true, variances are averaged over chromosomes
#'
#' @return if avg.over.chroms = TRUE, returns a data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Values "dX"
#' correspond to wavelet coefficients, and corresponding variances are associated
#'  with changes in localized average signal values over scale 2^(X-1). Values "sX" correspond to
#'  scaling coefficients, and corresponding variances are associated with localized averages
#'  on a scale of 2^X.  \cr
#'  \tab \cr
#'  \code{variance.X} \tab A weighted average across chromosomes of the magnitude of the
#'  wavelet variance of signal "X" at the specified level of the decomposition.
#'  Chromosomes without any coefficients at the specified level do not contribute to this average.
#'  For the level of the entire chromosome, the value is a weighted variance of chromosome signal
#'  means.
#'  \cr
#'  \tab \cr
#'  \code{propvar.X} \tab A weighted average across chromosomes of the proportion of the
#'  total variance of signal "X" contributed by the specified level of the decomposition.
#'  Values are then adjusted to account for the total proportion of genomic-wide variance contributed
#'  by differences in chromosome means. Thus, this has the interpretation of the total proportion
#'  of genome-wide variance of "X" contributed by the specified level.
#' @export
#'
#' @examples
#' # Simulated example - two white-noise signals "x" and "y" measured along two chromosomes.
#' # Chromosomes differ in their means by the same amount, and signal y has larger variance than x.
#' library(data.table)
#' x <- data.table(x = c(rep(0, 30), rep(10, 70)) + rnorm(100),
#'                 y = c(rep(0, 30), rep(10, 70)) + rnorm(100, 0, 10),
#'                 group = c(rep(1, 30), rep(2, 70))
#'                 position.id = c(1:30, 1:70))
#' plot(x[, c(x,y)])
#'
#' wv <- gnom_wav_var(x, chromosome = "group", signals = c("x","y"))
#' with(wv, plot(variance.y ~ level))
#' with(wv, plot(propvar.y ~ level))
#'
gnom_var_decomp <- function(x, chromosome, signals, rm.boundary=TRUE, avg.over.chroms = TRUE){
  m <- multi_modwts(data = x, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)

  # compute weights for each chromosome at each level of the decomposition
  # based on the number of wavelets present at that level on the given chromosome
  wav_weights <- m[, .(n.wavelets = .N), by = c(chromosome, "level")]

  # add NA weights for chromosomes missing certain levels
  for(g in unique(wav_weights$group)){
    absent <- setdiff(wav_weights[, level], wav_weights[group == g, level])
    wav_weights <- rbind(wav_weights,
                         list(group = g, level = absent, n.wavelets = rep(0,length(absent))))
  }

  wav_weights[, weight := n.wavelets/sum(n.wavelets), by = level][, n.wavelets := NULL]

  # compute wavelet variances using average of squared non-boundary wavelet coefficients
  wv <- m[, lapply(.SD, sample_var), .SDcols = paste0("coefficient.",signals), by = c(chromosome,"level")]
  varcols <- paste0("variance.", signals)
  setnames(wv, paste0("coefficient.",signals), varcols)

  # return chromosome results
  if (!avg.over.chroms) {
    return(wv)
  }

  # otherwise, return genome-wide results:

  totalvar <- x[, unlist(lapply(.SD, var)), .SDcols = signals]

  # ---- obtain chromosome-level variance
  # 1.  total signal means
  totalmeans <- x[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.",signals))

  # 2. chromosome means
  chrmeans <- x[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # 3. chromosome lengths as weights
  chrlen <- x[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)
  chrmeans[, names(totalmeans) := totalmeans]

  # 4. weighted variance
  chrvar <- vector(length=length(signals))
  names(chrvar) <- signals
  for(i in signals){
    j <- paste0("totalmean.",i)
    chrvar[i] <- chrmeans[, sum(weight*(get(i)- get(j))^2)]
  }

  # proportion of genomic variance from chromosome-level
  chrpropvar <- chrvar/totalvar

  # -----

  # merge with weights
  wv <- merge(wv, wav_weights, by = c(chromosome, "level"), all=T)

  # ----- calculate average magnitudes
  wv_mag <- wv[, lapply(.SD, weighted.mean, w = weight, na.rm=T),
               .SDcols = varcols, by = level]

  # ----- now calculate proportion of genomic variance***
  # 1. change NA values to zero to indicate that these scales don't contribute any variance for these chroms
  setnafill(wv, fill = 0, cols = varcols)

  # 2. calculate proportion of variance per chromosome
  propvar <- cbind(wv[, .(level)],
                   wv[, lapply(.SD, function(x){x/sum(x)}), .SDcols = varcols, by = chromosome])
  propvarcols <- paste0("propvar.",signals)
  setnames(propvar, varcols, propvarcols)

  # 3. average over chromosomes
  propvar <- propvar[, lapply(.SD, mean), .SDcols = propvarcols, by = level]

  allvardecomp <- merge(wv_mag, propvar)

  #*** we need to correct for the missing chromosome-level variance
  for(i in signals){
    chrvarcol <- paste0("propvar.", i)
    allvardecomp[, (chrvarcol) := get(chrvarcol)*(1-chrpropvar[i])]
  }

  # finally, add chromosome-level data to table

  chrtbl <- as.data.table(as.list(chrvar))
  setnames(chrtbl, signals, varcols)
  chrtbl <- cbind(chrtbl, as.data.table(as.list(chrpropvar)))
  setnames(chrtbl, signals, propvarcols)

  allvardecomp <- rbind(allvardecomp, data.table(level = chromosome, chrtbl))
  return(allvardecomp)
}
#
#
#
