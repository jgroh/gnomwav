#' Genome-wide variance decomposition
#'
#' Performs a scale-based decomposition of the variance in signal measurements along a genome
#' containing multiple chromosomes. Wavelet variances are obtained using the Maximal Overlap
#' Discrete Wavelet Transform with Haar wavelets. By default boundary coefficients are not
#' removed, but these can be ignored in the variance estimates. A weighted average
#' of wavelet variances is taken over chromosomes, and chromosome-level variance is also
#' calculated.
#'
#' @param data data.frame containing columns for chromosome id and signal measurements.
#' It is assumed that successive rows contain datapoints that are evenly spaced along the sequence.
#' If not, you will need to interpolate your signal to a grid of evenly spaced points.
#' @param chromosome character string, name of column containing chromosome id
#' @param signals character vector, names of columns containing the signals for which
#' variance decomposition is desired
#' @param rm.boundary logical, whether to remove boundary coefficients
#' @param avg.over.chroms logical, if FALSE, returns variance magnitudes per chromosome.
#' if TRUE, variance magnitudes are averaged over chromosomes,
#' along with 95% jacknife intervals from weighted block jacknife treating chromosomes as blocks.
#' Proportion of total genome-wide variance is also given.
#' Due to chromosomes having different lengths, error bars are not calculated for proportions of variance.
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
#'   \cr
#'   \tab \cr
#'  \code{variance.X.jack} \tab Weighted jacknife estimate of the variance. At each level of the decomposition,
#'  only chromosomes containing coefficients at that level are used in the jacknife.
#'   \cr
#'   \tab \cr
#'  \code{variance.X.jack.se} \tab Standard error of the jacknife estimate.
#' }
#' @export
#'
#' @examples
#' # Simulated example - two white-noise signals "x" and "y" measured along two chromosomes. Signal y has larger variance than x.
#'
#' library(data.table)
#' library(ggplot(2))
#' d <- data.table(chr = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000), y = rnorm(4000, sd=10))
#' gnom_var_decomp(d, chromosome = "chr", signals = c("x","y"))


# data <- data.table(group = c(rep(1,1000), rep(2,1000), rep(3,750), rep(4,500), rep(5,500), rep(6,250)), x = rnorm(4000))
# data[, y := x + rnorm(4000, mean=0, sd=2)]
# signals <- c('x', 'y'); chromosome <- 'group'; avg.over.chroms=T; rm.boundary <- F

gnom_var_decomp <- function(data, chromosome, signals, rm.boundary=FALSE, avg.over.chroms = TRUE){

  if(rm.boundary){
    warning("The scaling variance may be badly biased if boundary coefficients are removed. Consider comparing results to rm.boundary=F")
  }
  m <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)


  # ===== single chromosome
  if ( is.na(chromosome) ){
    wvd <- m[grepl("d", level), lapply(.SD, function(x){mean(x^2)}), .SDcols = paste0("coefficient.",signals), by = level]

    wvs <- m[grepl("s", level), lapply(.SD, function(x){mean(x^2) - mean(x)^2}), .SDcols = paste0("coefficient.",signals), by = level]
    wv <- rbind(wvd, wvs)

    varcols <- paste0("variance.", signals)
    setnames(wv, paste0("coefficient.",signals), varcols)
    return(wv)
  }

  # ===== multiple chromosomes

  n_chrom <- data[, length(unique(get(chromosome)))]
  wav_weights <- m[, .(n.wavelets = .N), by = c(chromosome, "level")]

  # number of classes of decomposition level
  n_class <- wav_weights[, .(nlevs = length(unique(level))), by = chromosome][, length(unique(nlevs))]

  if ( n_chrom > 1 && n_class > 1 ) {

    # add NA weights for chromosomes missing certain levels (only if the decomp levels actually differ across chromosomes)
    for(g in unique(wav_weights[, get(chromosome)])){
      absent <- setdiff(wav_weights[, level], wav_weights[get(chromosome) == g, level])

      abs_tbl <- data.table(chromosome = g,
                            level = absent,
                            n.wavelets = rep(0, length(absent) ) )
      setnames(abs_tbl,  "chromosome", chromosome)
      wav_weights <- rbind(wav_weights, abs_tbl)
    }
  }

  wav_weights[, weight := n.wavelets/sum(n.wavelets), by = level][, n.wavelets := NULL]

  # compute wavelet variances per chromosome using average of squared wavelet coefficients
  wvd <- m[grepl("d", level, fixed = T), lapply(.SD, function(x){mean(x^2)}), .SDcols = paste0("coefficient.",signals), by = c(chromosome,"level")]
  # scaling variances subtract the squared mean of coefficients
  wvs <- m[grepl("s", level, fixed = T), lapply(.SD, function(x){mean(x^2) - mean(x)^2}), .SDcols = paste0("coefficient.",signals), by = c(chromosome,"level")]
  wv <- rbind(wvd, wvs)

  varcols <- paste0("variance.", signals)
  setnames(wv, paste0("coefficient.",signals), varcols)

  # option to return chromosome results
  if (!avg.over.chroms) {
    return(wv)
  }

  # ----- otherwise, average across chroms to return genome-wide results:
  totalvar <- data[, unlist(lapply(.SD, var)), .SDcols = signals]

  # ---- obtain chromosome-level variance
  # 1.  total signal means
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.",signals))

  # 2. chromosome means
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # 3. chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)
  chrmeans[, names(totalmeans) := totalmeans]

  # 4. weighted variance
  chrvar <- vector(length=length(signals))
  names(chrvar) <- signals
  for(s in signals){
    k <- paste0("totalmean.",s)
    chrvar[s] <- chrmeans[, sum(weight*(get(s)- get(k))^2)]
  }

  # proportion of genomic variance from chromosome-level
  chrpropvar <- chrvar/totalvar

  # ----- Average wavelet variances across chromosomes

  # merge with weights
  wv <- merge(wv, wav_weights, by = c(chromosome, "level"), all=T) # all=T maintains levels in resulting dt not present on choromosomes

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

  # 3. average over chromosomes at each level, weight by chromosome length
  propvar <- merge(propvar, chrlen)
  propvar <- propvar[, lapply(.SD, weighted.mean, w=weight), .SDcols = propvarcols, by = level]

  allvardecomp <- merge(wv_mag, propvar)

  #*** we need to correct for the missing chromosome-level variance
  for(sig in signals){
    varcol <- paste0("propvar.", sig)
    allvardecomp[, (varcol) := get(varcol)*(1-chrpropvar[sig])]
  }

  # finally, add chromosome-level data to table

  chrtbl <- as.data.table(as.list(chrvar))
  setnames(chrtbl, signals, varcols)
  chrtbl <- cbind(chrtbl, as.data.table(as.list(chrpropvar)))
  setnames(chrtbl, signals, propvarcols)

  allvardecomp <- rbind(allvardecomp, data.table(level = chromosome, chrtbl))
  #return(allvardecomp)


  # ====== weighted jacknife =====
  # https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf
  # weighted jacknife of variance magnitudes to obtain confidence intervals.
  # Not straightforward to do confidence intervals of proportion of variance, due to differences in chromosome number across levels.
  # but confidence intervals of variance magnitudes can later be normalized to give a related interpretation.

  varcols_jack <- paste0(varcols, ".jack")
  varcols_jack_se <- paste0(varcols, ".jack.se")


  # compute separately for each level, as each level has potentially different number of chromosomes
  for(j in allvardecomp$level){
    if(j != chromosome){
      w_sub <- m[level == j]
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

    # don't do jacknife if too few chroms
    if(g < 4){
      for(vcol in varcols){
        allvardecomp[level==j, paste0(vcol,".jack") := NA][]
        allvardecomp[level==j, paste0(vcol,".jack.se") := NA][]
      }
    } else if (g >= 4){

      # leave-1-out estimates will go in here
      var_j <- data.table(jchr = jchrs, level = j)

      for(i in 1:length(jchrs)){
        # loop over chroms and leave each out in turn

        if(grepl('d',j,fixed=T)){
          # detail coefficient
         var_j[i, c(varcols)] <- w_sub[get(chromosome) != jchrs[i], lapply(.SD, function(x){mean(x^2)}), .SDcols = paste0("coefficient.",signals)]
        } else if(grepl('s',j,fixed=T)){
          # smooth coefficient
          var_j[i, c(varcols)] <- w_sub[get(chromosome) != jchrs[i], lapply(.SD, function(x){mean(x^2) - mean(x)^2}), .SDcols = paste0("coefficient.",signals)]
        } else if(j == chromosome){

          # ---- obtain chromosome-level variance
          # 1.  total signal means
          totalmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), .SDcols = signals]
          setnames(totalmeans_j, signals, paste0("totalmean.",signals))

          # 2. chromosome means
          chrmeans_j <- data[get(chromosome) != jchrs[i], lapply(.SD, mean), by = chromosome, .SDcols = signals]

          # 3. chromosome lengths as weights
          chrlen_j <- data[get(chromosome) != jchrs[i], .(weight = .N), by = chromosome]
          chrlen_j[, weight := weight/sum(weight)]
          chrmeans_j <- merge(chrmeans_j,chrlen_j)
          chrmeans_j[, names(totalmeans_j) := totalmeans_j]

          # 4. weighted variance
          for(s in signals){
            k <- paste0("totalmean.",s)
            var_j[i, paste0('variance.',s)] <- chrmeans_j[, sum(weight*(get(s)- get(k))^2)]
          }
        }
      }

      # compute jacknife estimates here
      h_j <- n/j_weights

      for(v in varcols){
        theta <- allvardecomp[level==j, get(v)]
        theta_j <- var_j[, get(v)]
        # jacknife estimate of the variance
        allvardecomp[level==j, paste0(v, ".jack") := g*theta - sum( ((n-j_weights)*theta_j)/n )][]

        # jacknife standard error of the estimator
        tau_j <- h_j*theta - (h_j - 1)*theta_j
        allvardecomp[level==j, paste0(v, ".jack.se") := sqrt((1/g)*sum((tau_j-theta_j)^2/(h_j-1)))][]
      }
    }
  }
  return(allvardecomp)
}
#
#
#
