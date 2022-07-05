#' Continuous Haar Wavelet
#'
#' @param y_H
#' @param j_H
#'
#' @return
#' @export
#'
#' @examples
haarCts <- function(y_H,j_H,k=1){

  # define wavelet support
  upper <- k*(2^j_H)
  lower <- upper - (2^j_H)
  mid <- upper - (2^j_H)/2

  return((y_H <= lower)*0 +
           (y_H > lower & y_H <= mid)*2^(-j_H/2) +
           (y_H > mid & y_H <= upper)*(-2^(-j_H/2)))
}


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
#' @param level
#' @param gen
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
wavelet_variance_equilbrium <- function(n.pop, n.sample, unit.scale, level, gen, alpha){
  genlist <- list()

  # loop over generations
  for(i in 1:length(gen)){
    t.gens <- as.numeric(gen[i])

    # make grid of parameters over which we evaluate the function
    grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, level=level, stringsAsFactors = F)
    grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size

    grd$gen <- rep(t.gens, nrow(grd)) # rep since we are inside the loop for a specific generation

    grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation

    for(q in 1:nrow(grd)){
      j <- grd[q,]$level
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



#' m
#' Transition probability matrix for haplotype state
#'
#' @param l_M
#' @param N_M
#' @param r_M
#' @param tau_M
#'
#' @return
#' @export
#'
#' @examples
trnstn_mat <- function(l_M,N_M,r_M=1/1024,tau_M){
  # l is a two vector of locus positions
  # N_M is population size
  # r_M is recombination distance per unit distance
  # tau_M is duration for which we want transition probabilities

  # only for l1 != l2
  # for l1==l2, haplotype state probs calculated separately in integrand funciton
  l <- as.numeric(l_M)
  L <- abs(l[2]-l[1])

  # matrix of right eigenvectors (A)
  A <- rbind(c(1,1),c(-1/(N_M*r_M*L),1))
  # matrix exponential of diagonal matrix of eigenvalues
  Dexp <- rbind( c( exp(-tau_M*((1 + N_M*r_M*L)/N_M)), 0), c(0,1))
  return(A %*% Dexp %*% solve(A))
}


#' wvBottleneckIntegrand
#'
#' @param x_I
#' @param J_I
#' @param Nvec_I
#' @param genvec_I
#' @param r_I
#' @param n.sample_I
#' @param alpha_I
#'
#' @return
#' @export
#'
#' @examples
wvBottleneckIntegrand <- function(x_I,J_I,Nvec_I,genvec_I,r_I=1/1024,n.sample_I,alpha_I=0.5){
  # x is a 2-vector of neutral locus positions

  # Nvec is number of haploid chromosomes, a vector through time
  # genvec are the time intervals, ordered toward the present. *These are durations of intervals.
  # Note these can be scalar values for a constant pop size.

  # r is recombination distance in M per unit distance of the signal
  # n.sample is the number of sampled chromosomes
  # alpha is the initial admixture proportion

  if(x_I[1] == x_I[2]){
    # l1==l2
    ii <- alpha_I

    # probability not coalesced:
    p_nc <- exp(-sum(genvec_I/Nvec_I))
    p_c <- 1 - p_nc
    ij <- alpha_I*p_c + (alpha_I^2)*p_nc

  } else {
    # l1! = l2: 2-state markov model for haplotype state probabilities

    # matrix product of transition probability matrices
    # leftmost matrix is starting from the present
    P <- m(l_M=x_I, N_M=Nvec_I[length(Nvec_I)], tau_M=genvec_I[length(genvec_I)])

    if(length(genvec_I) > 1){
      for(i in (length(genvec_I)-1):1){
        P <- P %*% trnstn_mat(l_M=x_I, N_M=Nvec_I[i], tau_M=genvec_I[i])
      }
    }

    ii_ij <- P %*% matrix(c(alpha_I,alpha_I^2), nrow=2)
    ii <- ii_ij[1,]
    ij <- ii_ij[2,]
  }

  return(haarCts(y_H=x_I[1],j_H=J_I)*haarCts(y_H=x_I[2],j_H=J_I)*(
    ((1/n.sample_I)*ii + ((n.sample_I-1)/n.sample_I)*ij) ) )
}

#' dblRiemann
#' helper function used for approximate integration
#'
#' @param scale_R
#' @param nMesh
#' @param Nvec_R
#' @param genvec_R
#' @param n.sample_R
#'
#' @return
#' @export
#'
#' @examples
dblRiemann <- function(scale_R, nMesh=100, Nvec_R, genvec_R, n.sample_R){
  # right Riemann sum approx of integral over a square region
  pnts <- seq(0,2^scale_R,length.out=nMesh+1)
  xy <-  expand.grid(pnts[-1], pnts[-1]) # take off 1st value for right riemann sum

  totalVol <- 0
  for(i in 1:nrow(xy)){
    x <- as.numeric(xy[i,])
    height <- wvBottleneckIntegrand(x, Nvec_I=Nvec_R, genvec_I=genvec_R, n.sample_I=n.sample_R, J_I=scale_R)
    totalVol <- totalVol + height*(2^scale_R/nMesh)^2
  }
  return(totalVol/2^scale_R)
}


#' wvBottleneck
#' operates on a vector of parameters in a table
#'
#' @param d
#' @param popSizeModel
#' @param epochs
#'
#' @return
#' @export
#'
#' @examples
wvBottleneck <- function(d, n.pop, epochs){
  epochCum <- c(0,cumsum(epochs))
  n.sample <- d[,n.sample]
  scl <- d[,level]
  g <- d[,gen]

  # construct Nvec and genvec based on gen
  t0 <- g - max(epochCum[epochCum < g]) # gives remaining duration of interval spent in first epoch toward past
  if(max(epochCum[epochCum < g]) == 0){
    Tau <- t0
    N <- n.pop[1]
  } else{
    tn <- epochs[1:(which(epochCum == max(epochCum[epochCum < g])) - 1)] # gets duration of subsequent intervals
    Tau <- c(tn, t0) # this vec
    N <- n.pop[1:(which(epochCum == max(epochCum[epochCum < g])))]
  }

  a <- dblRiemann(scale_R=scl, nMesh=100, Nvec_R=N, genvec_R=Tau, n.sample_R=n.sample)
  return(a)
}


wavelet_variance_general <- function(n.pop, epochs, n.sample, unit.scale, level, gen, alpha){

  grd1 <- data.table(expand.grid(gen=gen, n.sample=n.sample, level=level, stringsAsFactors = F))
  grd1[, variance := wvBottleneck(.SD, n.pop=n.pop, epochs=epochs), by = seq_len(nrow(grd1))][]
  return(grd1)
}




