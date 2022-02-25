cor_tbl <- function(data, chromosome, signals, rm.boundary = TRUE){
  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  cols <- paste0("coefficient.", signals)

  # wavelet correlations
  cor_tbl <- w[, .(cor = cor(get(cols[1]),get(cols[2]))), by =  level]

  # chromosome-scale correlation:
  totalmeans <- data[, lapply(.SD, mean), .SDcols = signals]
  setnames(totalmeans, signals, paste0("totalmean.", signals))
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = signals]

  # chromosome lengths as weights
  chrlen <- data[, .(weight = .N), by = chromosome]
  chrlen[, weight := weight/sum(weight)]
  chrmeans <- merge(chrmeans,chrlen)
  chrmeans[, names(totalmeans) := totalmeans]


  # 4. weighted chromosome-scale covariance
  totalmeancols <- paste0("totalmean.", signals)
  chrcov <- chrmeans[, mean(weight*(get(signals[1]) - get(totalmeancols[1]))*(get(signals[2]) - get(totalmeancols[2])))]

  wv <- gnom_var_decomp(data = data, chromosome = chromosome, signals = signals, avg.over.chroms = T)
  varcols <- paste0("variance.", signals)
  denom <- wv[level == chromosome, sqrt(get(varcols[1])*get(varcols[2]))]

  cor_tbl <- rbind(cor_tbl,
                  data.table(level = chromosome, cor = chrcov/denom))
  return(cor_tbl)

}



gnom_cor_decomp <- function(data, chromosome, signals, rm.boundary = TRUE){

  cor_n <- cor_tbl(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  setnames(cor_n, "cor", "cor_n")

  n_tot <- data[, length(unique(get(chromosome)))]

  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  nlev <- w[, .SD[, .(nlev = length(unique(get(chromosome))))], by = level]

  # this outputs a table of cors for all levels for each dropped chromosome
  cortbl_drop1 <- data[, cor_tbl(data[get(chromosome) != .BY],
                             chromosome=chromosome,
                             signals=signals,
                             rm.boundary = rm.boundary),
                     by = chromosome]
  cortbl_drop1 <- merge(cortbl_drop1, cor_n, by = "level")

  # jacknife bias-corrected estimates
  cortbl_drop1 <- merge(cortbl_drop1, nlev)
  cor_jack <- cortbl_drop1[, ps := nlev*(cor_n) - (nlev-1)*cor][, .(cor_jack = mean(ps)), by = level]

  # jacknife standard errors
  cor_se_jack <- cortbl_drop1[, .(se = sqrt(var(cor)/nlev)), by = level]

  output <- merge(cor_jack, cor_se_jack)
  output[, c("lower95ci", "upper95ci") := .(cor_jack - 1.96*se, cor_jack + 1.96*se)][, se := NULL][]
  output <- merge(cor_n, output, by = "level")

  return(output)
}


jacknife_rsqrd <- function(data, yvar, xvars, chromosome){

  f <- as.formula(paste(yvar, paste(xvars, collapse=" + "), sep=" ~ "))

  # model using all data
  z <- summary(lm(f, data = data))

  # estimate of r squared
  theta_n <- z$r.squared

  n <- length(data[, unique(get(chromosome))])

  # don't construct jacknife CIs if fewer than X # of chromosome
  if ( n < 3){
    rsqrd <- theta_n
    lower <- as.double(NA)
    upper <- as.double(NA)

  } else { # construct CIs

    # jacknife estimates
    theta_i <- vector()
    ps <- vector()
    cnt <- 0

    # note: include some check here to make sure there are sufficiently many chromosomes
    for(i in data[, unique(get(chromosome))]){
      cnt <- cnt + 1
      theta_i <- summary(lm(f, data = data[get(chromosome) != i,]))$r.squared
      ps[cnt] <- n*theta_n - (n-1)*theta_i
    }

    rsqrd <- mean(ps) # bias-corrected point estimate

    # 95% confidence intervals
    se <- sqrt(var(ps)/n)
    lower <- rsqrd - 1.96*se
    upper <- rsqrd + 1.96*se

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
