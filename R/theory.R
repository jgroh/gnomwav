# Expected wavelet variance: approximate integrand ------------------------------------------
# assumes infinite population
#' Title
#'
#' @param x
#' @param u
#' @param n.sample
#' @param alpha
#' @param t.gens
#'
#' @return
#'
#' @examples
wav_var_approx <- function(x, u, n.sample, alpha, t.gens) {
  (1/n.sample)*(
    alpha*exp(-t.gens*u*abs(x[2]-x[1])) +
      (alpha^2)*(1-exp(-t.gens*u*abs(x[2]-x[1])))
  ) + ((n.sample-1)/n.sample)*alpha^2
}

# Expected wavelet variance: exact  integrand ------------------------------------------

#' Title
#'
#' @param x
#' @param expected.crossovers.per.unit.dist
#' @param n.pop
#' @param n.sample
#' @param alpha
#' @param t.gens
#'
#' @return
#'
#'
#' @examples
wav_var_exact <- function(x, expected.crossovers.per.unit.dist,
                          n.pop, n.sample, alpha, t.gens) {
  u <- expected.crossovers.per.unit.dist

  v <- n.pop*u*abs(x[2]-x[1])
  w <- exp(-(t.gens/n.pop)*(1+v))
  (
    (1/n.sample)*(
      alpha*(1+v*w)/(1+v) + alpha^2*(v*(1-w))/(1+v)
    )
    +
      ((n.sample-1)/n.sample)*(
        alpha*(1-w)/(1+v) + alpha^2*(v+w)/(1+v))
  )
}

#' Title
#'
#' @param n.pop
#' @param n.sample
#' @param scale
#' @param gen
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
wavelet_variance_equilbrium <- function(n.pop, n.sample, unit.scale, scale, gen, alpha){
  genlist <- list()

  # loop over generations
  for(i in 1:length(gen)){
    t.gens <- as.numeric(gen[i])

    # make grid of parameters over which we evaluate the function
    grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, scale=scale, stringsAsFactors = F)
    grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size

    grd$gen <- rep(t.gens, nrow(grd)) # rep since we are inside the loop for a specific generation

    grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation

    for(q in 1:nrow(grd)){
      j <- grd[q,]$scale
      ns <- grd[q,]$n.sample
      np <- grd[q,]$n.pop

      if(n.pop == Inf){ # use infinite population approximation
        part1 <- adaptIntegrate(wav_var_approx, n.sample = ns, expected.crossovers.per.unit.dist=unit.scale, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,0),
                                upperLimit = c(2^(j-1),2^(j-1)))
        part2 <- adaptIntegrate(wav_var_approx, n.sample = ns, expected.crossovers.per.unit.dist=unit.scale, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,2^(j-1)),
                                upperLimit = c(2^(j-1),(2^j)))
      } else { # use exact formula
        part1 <- adaptIntegrate(wav_var_exact, n.sample = ns, n.pop = np, expected.crossovers.per.unit.dist=unit.scale, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,0),
                                upperLimit = c(2^(j-1),2^(j-1)))
        part2 <- adaptIntegrate(wav_var_exact, n.sample = ns, n.pop = np, expected.crossovers.per.unit.dist=unit.scale, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,2^(j-1)),
                                upperLimit = c(2^(j-1),(2^j)))
      }
      grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
    }

    genlist[[i]] <- grd
  }

  return(data.table(do.call(rbind.data.frame, genlist)))
}










