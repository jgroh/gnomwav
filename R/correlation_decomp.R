gnom_cor_decomp <- function(data, chromosome, signals, rm.boundary = TRUE){

  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)

  # compute correlations
  cols <- paste0("coefficient.", signals)
  cor_tbl <- w[, .(cor = cor(get(cols[1]),get(cols[2]))), by =  level]

  # compute chromosome-scale correlation
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.",signals))

  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)
  chrmeans[, names(totalmeans) := totalmeans]

  # 4. weighted covariance
  totalmeancols <- paste0("totalmean.", signals)
  chrcov <- chrmeans[, mean(weight*(get(signals[1]) - get(totalmeancols[1]))*(get(signals[2]) - get(totalmeancols[2])))]

  wv <- gnom_var_decomp(data = data, chromosome = chromosome, signals = signals, avg.over.chroms = T)
  varcols <- paste0("variance.", signals)
  denom <- wv[level == chromosome, sqrt(get(varcols[1])*get(varcols[2]))]

  return(rbind(cor_tbl,
               data.table(level = chromosome, cor = chrcov/denom)))

}


jacknife_rsqrd <- function(data, yvar, xvars, chromosome){

  f <- as.formula(paste(yvar, paste(xvars, collapse=" + "), sep=" ~ "))

  # model using all data
  z <- summary(lm(f, data = data))
  theta_n <- z$adj.r.squared

  # estimate of r squared
  theta_n <- z$adj.r.squared

  # don't construct jacknife CIs if fewer than X # of chromosome
  if ( data[, length(unique(get(chromosome)))] < 3){
    rsqrd <- theta_n
    lower <- as.double(NA)
    upper <- as.double(NA)

  } else { # construct CIs

    # pseudovalues
    ps <- vector()
    cnt <- 0
    n <- length(data[, unique(get(chromosome))])

    # note: include some check here to make sure there are sufficiently many chromosomes
    for(i in data[, unique(get(chromosome))]){
      cnt <- cnt + 1
      theta_i <- summary(lm(f, data = data[get(chromosome) != i,]))$adj.r.squared
      ps[cnt] <- n*theta_n - (n-1)*theta_i
    }

    # 95% confidence interval
    #rsqrd <- mean(ps)
    rsqrd <- theta_n
    ps_var <- var(ps)

    lower <- rsqrd - 1.96*sqrt(ps_var/n)
    upper <- rsqrd + 1.96*sqrt(ps_var/n)

  }

  return( c(list(rsqrd = rsqrd, ci95_lower = lower, ci95_upper = upper)))

}

wvlt_lm_rsqrd <- function(data, yvar, xvars, chromosome, rm.boundary = TRUE){

  w <- multi_modwts(data = data, chromosome = chromosome, signals = c(xvars,yvar), rm.boundary = rm.boundary)

  # calculate r squared and 95% confidence intervals across scales
  rsqrd_tbl <- w[, jacknife_rsqrd(.SD,
                                  chromosome = chromosome,
                                  xvars = paste0("coefficient.",xvars),
                                  yvar = paste0("coefficient.",yvar)),
                 by = level]

  # calculate r squared and 95% for chromosome-level
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = c(xvars,yvar)]

  rsqrd_chr <- cbind(level = chromosome,
                           chrmeans[, jacknife_rsqrd(.SD,
                                                     chromosome = chromosome,
                                                     xvars = xvars,
                                                     yvar = yvar)])
  return(rbind(rsqrd_tbl, rsqrd_chr))
}

#load("~/workspace/selection-against-introgression/datasets/swordtail_TOTO_ACUA/ACUA_2018/gnomP.RData")
#gnomP <- gnomP[ID==ID[1]]


