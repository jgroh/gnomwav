# helper function to return vector of wavelet levels as named by waveslim
wt_levels <- function(data, chromosome){

  # define maximum level of decomposition based on longest chromosome
  maxLev <- data[, .(n_obs = .N), by = chromosome][, max(floor(log2(n_obs)))]

  # names expected from output of waveslim::modwt
  return(c(paste0("d", 1:maxLev), paste0("s", maxLev)))
}

#' modwt of a signal on one chromosome of a genome with multiple chromosomes
#'
#' Helper function for dijoint_modwts below.
#' Wrapper for modwt function from waveslim. Sets coefficients to NA for levels
#' that are not present on the focal chromosome but that are present on others.
#'
#' @param u a vector of signal measurements
#'
#' @return a data.table with 2 columns, level of the decomposition,
#' coefficient values
#'
#' @import data.table
#'
#'
multi_modwt_1var <- function(u, all_cols, rm.boundary = F){
  a <- waveslim::modwt(x = u,wf = "haar",n.levels = floor(log2(length(u))))

  if(rm.boundary){
    a <- waveslim::brick.wall(wf = "haar", x = a)
  }

  b <- suppressWarnings(melt(as.data.table(a[names(a)]),
                          variable.name = "level",
                          value.name = "coefficient"))

  # I think no longer needed:
  #b[setdiff(all_cols, names(b))] <- as.double(NA) # coefficients for higher levels not present set to NA

  return(b)
}


#' Maximal Overlap Discrete Wavelet Transform of signals measured along
#' a set of disjoint sequences.
#'
#' Wrapper function for functions from package `waveslim`.
#' Returns a table of coefficients from the Maximal Overlap Discrete Wavelet Transform (MODWT)
#' performed on multiple signal measurements (i.e. ancestry state and recombination rate)
#' along multiple disjoint sequences (i.e. chromosomes). By default the 'brick wall' boundary coefficient
#' is applied, meaning that boundary coefficients are removed from the output.
#'
#' @param data a data.frame or data.table containing a column for chromosome ID and column(s) with signal values.
#' @param chromosome character string, the name of the column with chromosome ID
#' @param signals character vector, names of the columns with signal values
#'
#' @return data.table containing:\tabular{ll}{
#' \code{level} \tab The level of the wavelet decomposition. Levels "d?" correspond to wavelet coefficients.
#'  Levels "s?" correspond to scaling coefficients. \cr
#'  \tab \cr
#'  \code{coefficient.X} \tab Coefficient values for the signal named "X". \cr
#'  \tab \cr
#'  \code{position.id} \tab  End position of the wavelet or scaling filter that produced this coefficient.
#'  As boundary coefficients are removed, the first position.id is the end position of the wavelet that starts
#'  at the beginning of the sequence.  \cr
#' }
#'
#' @import data.table
#' @export
#'
#' @examples tbd
#'
multi_modwts <- function(data, chromosome, signals, rm.boundary=F){
  setDT(data)
  d <- copy(data)

  if (is.na(chromosome)){
    d2 <- d[, multi_modwt_1var(u = get(signals[1]), all_cols = allcols, rm.boundary = rm.boundary)]
    d2[, position.id := seq_len(.N), by = "level"]

    setnames(d2, "coefficient", paste0("coefficient.", signals[1]))

    if (length(signals) > 1) {
      for(i in 2:length(signals)){
        temp <- d[, multi_modwt_1var(get(signals[i]), all_cols = allcols, rm.boundary = rm.boundary)]
        temp[, position.id := seq_len(.N), by = "level"]
        setnames(temp, "coefficient", paste0("coefficient.", signals[i]))
        d2 <- merge(d2, temp, by = c("level", "position.id"))
      }
    }
  } else{

    allcols <- wt_levels(d, chromosome)

    # modwt on first signal
    d2 <- d[, multi_modwt_1var(u = get(signals[1]), all_cols = allcols, rm.boundary = rm.boundary), by = chromosome]
    d2[, position.id := seq_len(.N), by = c(chromosome, "level")]
    setnames(d2, "coefficient", paste0("coefficient.", signals[1]))

    # modwt for other signals and combine result
    if (length(signals) > 1) {
      for(i in 2:length(signals)){
        temp <- d[, multi_modwt_1var(get(signals[i]), all_cols = allcols, rm.boundary = rm.boundary), by = chromosome]
        temp[, position.id := seq_len(.N), by = c(chromosome, "level")]
        setnames(temp, "coefficient", paste0("coefficient.", signals[i]))
        d2 <- merge(d2, temp, by = c(chromosome, "level", "position.id"))
      }
    }
  }


  return(na.omit(d2, cols = paste0("coefficient.", signals)))
}

