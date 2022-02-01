gnom_cor_decomp <- function(data, chromosome, signals, rm.boundary = TRUE){
  w <- multi_modwts(data = data, chromosome = chromosome, signals = signals, rm.boundary = rm.boundary)
  wav_weights <- w[, .(n.wavelets = .N), by = c(chromosome, "level")]
  wav_weights[, weight := n.wavelets/sum(n.wavelets), by = level]

  wv <- gnom_var_decomp(data = data, chromosome = chromosome, signals = signals, avg.over.chroms = T, rm.boundary = rm.boundary )
  cols <- paste0("coefficient.", signals)

  cov_tbl <- w[, .(cov = mean(get(cols[1])*get(cols[2]))), by =  level]

  varcols <- paste0("variance.",signals)
  cov_tbl <- merge(cov_tbl, wv, all=T)

  corcol <- paste0("pearson.cor.", paste(signals, collapse ="."))
  cor_tbl <- cov_tbl[, .(level, cor = cov/sqrt(get(varcols[1])*get(varcols[1])))]
  setnames(cor_tbl, "cor", corcol)

  # chromosome scale
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
  chrcov <- chrmeans[, mean((get(signals[1]) - get(totalmeancols[1]))*(get(signals[2]) - get(totalmeancols[2])))]

  cor_tbl[level == chromosome, (corcol) := chrcov][]

  return(cor_tbl)

}

jacknife_rsqrd <- function(data, yvar, xvars, chromosome){

  f <- as.formula(paste(yvar, paste(xvars, collapse=" + "), sep=" ~ "))

  # estimate with all data
  theta_n <- summary(lm(f, data = data))$adj.r.squared

  # pseudovalues
  ps <- vector()
  cnt <- 0
  n <- length(data[, unique(get(chromosome))])

  for(i in data[, unique(get(chromosome))]){
    cnt <- cnt + 1
    theta_i <- summary(lm(f, data = data[get(chromosome) != i,]))$adj.r.squared
    ps[cnt] <- n*theta_n - (n-1)*theta_i
  }

  # 95% confidence interval
  ps_mean <- mean(ps)
  ps_var <- var(ps)

  low <- ps_mean - 1.96*sqrt(ps_var/n)
  up <- ps_mean + 1.96*sqrt(ps_var/n)

  return(list(rsqrd = ps_mean, ci95_lower = low, ci95_upper = up))

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



