
#data <- data.table(x = rnorm(100), y = rnorm(100),
#                   z = 1:100, chromosome = c(rep(1,70), rep(2,30)))


brickWallModwt <- function(data, chromosome, signals){
  setDT(data)
  d <- copy(data)
  setnames(d, get(chromosome), "group")

  # define maximum level of decomposition based on longest chromosome
  maxLev <- d[, .(n_obs = nrow(.SD)), by = group][, max(floor(log2(n_obs)))]

  # names expected from output of waveslim::modwt
  allCols <- c(paste0("d", 1:maxLev), paste0("s", maxLev))

  #
  modwtOneVar <- function(u){
    b <- waveslim::brick.wall(wf = "haar",
                              x = waveslim::modwt(x = u,
                                                  wf = "haar",
                                                  n.levels = floor(log2(length(u)))))
    b[setdiff(allCols, names(b))] <- as.double(NA) # coefficients for higher levels not present set to NA

    b <- suppressWarnings(melt(as.data.table(b[allCols]),
                               variable.name = "level",
                               value.name = "coefficient"))
    b[, position_id := seq_len(.N)]
  }


  d2 <- d[, modwtOneVar(get(signals[1])), by = group]
  setnames(d2, "coefficient", paste0("coefficient.", signals[1]))

  for(i in 2:length(signals)){
    temp <- d[, f(get(signals[i])), by = group]
    setnames(temp, "coefficient", paste0("coefficient.", signals[i]))
    d2 <- merge(d2, temp)
  }

  return(d2)
}

