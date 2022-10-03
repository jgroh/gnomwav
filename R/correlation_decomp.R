
cov_tbl <- function(data, chromosome, signals, rm.boundary = TRUE){
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



cor_tbl <- function(data, chromosome, signals, rm.boundary = TRUE){

  # get modwt coefficients
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  cols <- paste0("coefficient.", signals)

  # wavelet 'correlations' these are not quite correlations bc we don't subtract off the product of the means

  cor_tbl <- w[, .(cor = mean(get(cols[1])*get(cols[2]))/ (sqrt(mean(get(cols[1])^2)*mean(get(cols[2])^2)))), by = level]

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
                  data.table(level = chromosome, cor = chrcor))
  return(cor_tbl)
}



gnom_cor_decomp <- function(data, chromosome, signals, rm.boundary = TRUE){

  cor_n <- cor_tbl(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  setnames(cor_n, "cor", "cor_n")

  n_tot <- data[, length(unique(get(chromosome)))]

  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  nlev <- w[, .SD[, .(nlev = length(unique(get(chromosome))))], by = level]
  nlev <- rbind(nlev, data.table(level = chromosome, nlev = n_tot))

  # this outputs a table of cors for all levels for each dropped chromosome
  cortbl_drop1 <- data[, cor_tbl(data[get(chromosome) != .BY],
                             chromosome=chromosome,
                             signals=signals,
                             rm.boundary = rm.boundary),
                     by = chromosome]
  cortbl_drop1 <- merge(cortbl_drop1, cor_n, by = "level")

  # jacknife bias-corrected estimates
  cortbl_drop1 <- merge(cortbl_drop1, nlev, all = T)

  cor_jack <- cortbl_drop1[, ps := nlev*(cor_n) - (nlev-1)*cor][, .(cor_jack = mean(ps)), by = level]

  # jacknife standard errors
  cor_se_jack <- cortbl_drop1[, .(se = ifelse(nlev[1] != 1, sqrt(var(cor)/nlev[1]), as.double(NA))), by = level]

  output <- merge(cor_jack, cor_se_jack)
  output[, c("lower95ci", "upper95ci") := .(cor_jack - 1.96*se, cor_jack + 1.96*se)][, se := NULL][]
  output <- merge(cor_n, output, by = "level")

  return(output)
}



jacknife_lm <- function(dt, y, x, chromosome){

  f <- as.formula(paste(y, paste(x, collapse=" + "), sep=" ~ "))

  # model using all data
  z <- summary(lm(f, data = dt))

  # linear model coeffs
  lm_coeffs <- z$coefficients[-1,1]

  # estimate of r squared
  rsqrd_n <- z$r.squared
  names(rsqrd_n) <- "rsqrd"

  estimates <- unlist(c(lm_coeffs, rsqrd_n) )

  n <- length(dt[, unique(get(chromosome))])

  # don't construct jacknife CIs if fewer than X # of chromosome
  if ( n < 4){
    nms <- names(estimates)
    output <- data.table(estimates)
    output[, variable := as.character(nms)]
    output[, jn_se := as.double(NA)]
    output[, jn_bc_estimate := as.double(NA)][]
    setnames(output, "estimates", "direct_estimate")
    setcolorder(output, c("variable", "direct_estimate", "jn_bc_estimate", "jn_se"))

  } else { # construct jacknife CIs for rsquared and slopes

    jack_estimates <- data.table()

    # note: include some check here to make sure there are sufficiently many chromosomes
    for(i in dt[, unique(get(chromosome))]){

      z <- summary(lm(f, data = dt[get(chromosome) != i,]))

      jack_estimates <- rbind(jack_estimates, c(as.list(z$coefficients[-1, 1]),
                                      rsqrd = z$r.squared),
                              fill=T)
    }


    pseudovals <- jack_estimates[, lapply(seq_along(names(.SD)),
                            function(i, y, nm){ n*estimates[nm[i]] - (n-1)*y[[i]]  },
                            y = .SD, nm = names(.SD)),
                   .SDcols = c(x, "rsqrd")]
    setnames(pseudovals, new = c(x, "rsqrd"))

    # bias-corrected estimates: mean of pseudovalues
    jack_bc_estimates <- pseudovals[, lapply(.SD, mean)]

    output <- melt(jack_bc_estimates, measure.vars = names(jack_bc_estimates), value.name = "jn_bc_estimate")
    output[, variable := as.character(variable)]
    # jacknife standard errors
    jack_se <- pseudovals[, lapply(.SD, function(x){ sqrt(var(x)/n) } )]
    output <- merge(output, melt(jack_se, measure.vars = names(jack_se), value.name = "jn_se"))

    output[, direct_estimate := estimates][]
    setcolorder(output, c("variable", "direct_estimate", "jn_bc_estimate", "jn_se"))

  }

  return(output)

}


#data <- data.table(x = rnorm(100), y = rnorm(100), z = rnorm(100), chr = rep(c(1,2,3,4,5), 20))
#setkey(data, chr)
#xvars <- c("x", "y")
#yvar <- "z"
#chromosome <- "chr"
#dt <- data
#rm.boundary <- T

wvlt_lm <- function(data, yvar, xvars, chromosome, rm.boundary = TRUE){

  w <- multi_modwts(data = data, chromosome = chromosome, signals = c(xvars,yvar), rm.boundary = rm.boundary)

  # calculate r squared and 95% confidence intervals across scales
  rsqrd_tbl <- w[, jacknife_lm(dt = .SD,
                                  chromosome = chromosome,
                                  x = paste0("coefficient.",xvars),
                                  y = paste0("coefficient.",yvar)),
                 by = level]
  rsqrd_tbl[, variable := gsub("coefficient.","", variable)]

  # calculate r squared and 95% for chromosome-level
  chrmeans <- data[, lapply(.SD, mean), by = chromosome, .SDcols = c(xvars,yvar)]

  rsqrd_chr <- cbind(level = chromosome,
                           chrmeans[, jacknife_lm(.SD,
                                                     chromosome = chromosome,
                                                     x = xvars,
                                                     y = yvar)])
  return(rbind(rsqrd_tbl, rsqrd_chr))
}
