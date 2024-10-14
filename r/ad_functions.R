## ################# FUNCTION DEFINITIONS

build_ad <- function(data, model_id = 4, data_id = 6) {

  # libraries
  library(mgcv)
  library(rstan)

  if (!model_id %in% c(1:10)) stop("model_id incorrectly specified")

  if (model_id == 1) model.txt <- "ptens_quar_nourep.stan"
  if (model_id == 2) model.txt <- "ptens_quar_nourep_aidiag1.stan"
  if (model_id == 3) model.txt <- "ptens_quar_agediag4_nourep.stan"
  if (model_id == 4) model.txt <- "ptens_quar_agediag5_nourep.stan"
  if (model_id == 5) model.txt <- "ptens_quar_agediag_nourep.stan"
  if (model_id == 6) model.txt <- "ptens_quar_agediag1_nourep.stan"
  if (model_id == 7) model.txt <- "ptens_quar_agediag2_nourep.stan"
  if (model_id == 8) model.txt <- "ptens_quar_agediag3_nourep.stan"
  if (model_id == 9) model.txt <- "ptens_quar_agediag1_rita.stan"

  spl.file.new <- here("bug/ptens_10_8.bug")

  ## mgcv splines specifications
  ## dummy data
  tmp <- runif(data$nquar * data$nage)

  yrs <- 1:data$nquar ## yrs is actually quarters, but convenient to call it like his
  ages <- 1:data$nage

  ## Creating a dataset for jagam
  jags.data <- list(
    age = rep(ages, each = data$nquar),
    yrs = rep(yrs, times = data$nage),
    D = tmp
  )

  m.list <- list(c(2, 1), c(2, 1))
  jagam.out <- jagam(D ~ te(yrs, age, bs = c("ps", "ps"), k = c(10, 8), m = m.list),
    family = gaussian, data = jags.data,
    file = spl.file.new, diagonalize = TRUE
  )

  # append to the data
  data$X <- jagam.out$jags.data$X
  data$ninfpars <- ncol(jagam.out$jags.data$X)
  data$S1 <- jagam.out$jags.data$S1

  if (model_id %in% c(2:6)) pars_save <- c("beta", "lambda", "vardelta", "d", "delta1", "delta2", "delta3", "delta4", "alpha", "lat_arr") ## saving age-specific intercept dx parameters
  if (model_id %in% c(8)) pars_save <- c("beta", "lambda", "lambda_d1", "lambda_d2", "lambda_d3", "lambda_d4", "d", "delta1", "delta2", "delta3", "delta4")
  if (model_id ==9) pars_save <- c("beta", "lambda", "vardelta", "d", "delta1", "delta2", "delta3", "delta4", "delta5", "alpha", "lat_arr") ## saving age-specific intercept dx parameters

  stan_model <- stan_model(file = here("stan", model.txt))

  return(list(
    "stan_data" = data,
    "stan_model" = stan_model,
    "pars_save" = pars_save
  ))
}

### Function to produce consistent output for each scenario
summary.fct <- function(x) {
  ## gives 95% credible intervals
  ## median and mode
  out <- rep(NA, 4)
  out[1] <- quantile(x, probs = 0.025)
  out[2] <- quantile(x, probs = 0.5)
  out[3] <- mean(x)
  out[4] <- quantile(x, probs = 0.975)
  names(out) <- c("2.5%", "50%", "mean", "97.5%")
  out
}

### Function to get posterior functions of incidence
### i.e incidence by age group
### from log-infections (per each posterior iteration)
### This function can aggregates incidence by year, depending on qt.flag
infs.fct <- function(infs, age.st, age.end, nt, type = "time", qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (!type %in% c("time", "age")) stop("Type can only be 'time' or 'age'")
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (qt.flag) qt.indx <- rep(1:(nt / 4), each = 4)
  if (type == "time") {
    if (!qt.flag) {
      call <- "rowSums(matrix(x, nr=nt)[,a1:a2])"
    } else {
      call <- "tapply(rowSums(matrix(x, nr=nt)[,a1:a2]), qt.indx, sum)"
    }
  } else {
    call <- "colSums(matrix(x, nr=nt)[,a1:a2])"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  apply(x, 2, summary.fct)
}

### Function to get posterior functions of incidence
### i.e incidence by age group
### from log-infections (per each posterior iteration)
### This plots always aggregated QUARTERLY incidence curve
infs.fct.qt <- function(infs, age.st, age.end, nt, qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    infs <- infs[, rep(1:ncol(infs), each = 4)] / 4 ## making yearly infs constant in quarters
    nt <- nt * 4 ## nt is quarter also for yearly model
  }
  call <- "rowSums(matrix(x, nr=nt)[,a1:a2])"
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  apply(x, 2, summary.fct)
}

### Incidence function... on the basis of the denominator
### provided by Anne... I do not have a distribution so dividing by median
### (this somehow underestimates uncertainty)
inc.fct <- function(infs, age.cl, last.yr, nt, qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (!age.cl %in% c("Tot", "15-34", "35-44", "45+")) stop("age.cl can only be 'Tot','15-34', '35-44','45+' ")
  if (age.cl == "Tot") {
    a1 <- 1
    a2 <- 52
    den <- N.msm.tot
  }
  if (age.cl == "15-34") {
    a1 <- 1
    a2 <- 20
    den <- N.msm.1534
  }
  if (age.cl == "35-44") {
    a1 <- 21
    a2 <- 30
    den <- N.msm.3544
  }
  if (age.cl == "45+") {
    a1 <- 31
    a2 <- 52
    den <- N.msm.45
  }

  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2]) / den "
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum) / den"
  }

  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  apply(x, 2, summary.fct)
}

### Peak function
### Function aimed at identifying the peak in infection in the last x years
peak.fct <- function(infs, age.st, age.end, nt, last.yr = 10, qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (last.yr > (nt / 4)) stop("Specified last year is too large")
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2])"
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum)"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  ### yields maximum in the last years
  apply(x, 1, which.max)
}

### Function for the shape of incidence curve
## in the last two years only increasing and decreasing
## in other years...
## monotonically increasing, inc (not monotone), dec (not monotone), monotonic decrease
shape.fct <- function(x) {
  ## x arbitrary vector...
  ## do not consider flat case here as not really the case...
  l <- length(x)
  o <- NA
  if (l <= 1) {
    stop("length of x must be greater than 1")
  } else if (l == 2) {
    if (x[1] < x[l]) o <- "inc" else o <- "dec"
  } else {
    if (x[1] < x[l]) {
      if (all(x == cummax(x))) {
        o <- "mon inc"
      } else {
        o <- "inc (not mon)"
      }
    } else {
      if (all(x == cummin(x))) {
        o <- "mon dec"
      } else {
        o <- "dec (not mon)"
      }
    }
  }
  return(o)
}

trend.fct <- function(infs, age.st, age.end, nt, last.yr = 3, qt.flag = FALSE) {
  if (last.yr > (nt / 4)) stop("Specified last year is too large")
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2])"
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum)"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  apply(x, 1, shape.fct)
}


### Function to get posterior functions of goodness of fit
### i.e HIV/AIDS/CD4 by age group
### from HIV/AIDS/CD4 (per each posterior iteration)
gof.fct <- function(type = "HIV", age.st, age.end, HIV, AIDS, CD4.st1, CD4.st2, CD4.st3, CD4.st4, nt, qt.flag) { ## infs is a matrix with log-infections for each posterior sample in each row
  if (!type %in% c("HIV", "AIDS", "CD4.st1", "CD4.st2", "CD4.st3", "CD4.st4")) stop("Type can only be 'HIV','AIDS','CD4.st1','CD4.st2','CD4.st3','CD4.st4'")
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (qt.flag) qt.indx <- rep(1:(nt / 4), each = 4)
  call <- paste(type, "[,a1:a2,]", sep = "")
  y <- t(apply(eval(parse(text = call)), 3, rowSums))
  if (qt.flag) y <- t(apply(y, 1, function(x) tapply(x, qt.indx, sum)))
  apply(y, 2, summary.fct)
}

### Function to deal with EXPECTED prevalence
prev.fct <- function(prev.arr, state, age.st, age.end) { ## infs is a matrix with log-infections for each posterior sample in each row
  if (!state %in% c("all", 1:4)) stop("State can only be 'all' or 1:4")
  a1 <- age.st - 14
  a2 <- age.end - 14

  if (state == "all") {
    prev <- prev.arr[, , 1, ] + prev.arr[, , 2, ] + prev.arr[, , 3, ] + prev.arr[, , 4, ]
  } else {
    prev <- prev.arr[, , state, ]
  }

  x <- t(apply(prev[, a1:a2, ], 3, rowSums))
  apply(x, 2, summary.fct)
}

## Function to get output in right format
## For C++ functions
d.mat.iter <- function(fit, iter, age.dx.flag, age.dx.spl.flag, rita.flag, nt) {
  if (!age.dx.flag) {
    d.mat <- matrix(NA, nr = 4, nc = nt)
    d.mat[1, ] <- fit[iter, d1.ind]
    d.mat[2, ] <- fit[iter, d2.ind]
    d.mat[3, ] <- fit[iter, d3.ind]
    d.mat[4, ] <- fit[iter, d4.ind]
  } else {
    if (!age.dx.spl.flag) {
     if (rita.flag) {
      d.mat <- array(NA, dim = c(5, nt, 4))
      d.mat[5, , 1] <- fit[iter, d5.ind.y]
      d.mat[5, , 2] <- fit[iter, d5.ind.my]
      d.mat[5, , 3] <- fit[iter, d5.ind.mo]
      d.mat[5, , 4] <- fit[iter, d5.ind.o]
    } else {
      d.mat <- array(NA, dim = c(4, nt, 4)) ## array to be consistent with gof_iter_fct in stancppfcts_simpleperv_agediag.cpp
      }
      d.mat[1, , 1] <- fit[iter, d1.ind.y]
      d.mat[2, , 1] <- fit[iter, d2.ind.y]
      d.mat[3, , 1] <- fit[iter, d3.ind.y]
      d.mat[4, , 1] <- fit[iter, d4.ind.y]
      d.mat[1, , 2] <- fit[iter, d1.ind.my]
      d.mat[2, , 2] <- fit[iter, d2.ind.my]
      d.mat[3, , 2] <- fit[iter, d3.ind.my]
      d.mat[4, , 2] <- fit[iter, d4.ind.my]
      d.mat[1, , 3] <- fit[iter, d1.ind.mo]
      d.mat[2, , 3] <- fit[iter, d2.ind.mo]
      d.mat[3, , 3] <- fit[iter, d3.ind.mo]
      d.mat[4, , 3] <- fit[iter, d4.ind.mo]
      d.mat[1, , 4] <- fit[iter, d1.ind.o]
      d.mat[2, , 4] <- fit[iter, d2.ind.o]
      d.mat[3, , 4] <- fit[iter, d3.ind.o]
      d.mat[4, , 4] <- fit[iter, d4.ind.o]
    } else {
      d.mat <- array(NA, dim = c(4, nt, nage))
      d.mat[1, , ] <- matrix(fit.mat[iter, d1.ind], nr = nt)
      d.mat[2, , ] <- matrix(fit.mat[iter, d2.ind], nr = nt)
      d.mat[3, , ] <- matrix(fit.mat[iter, d3.ind], nr = nt)
      d.mat[4, , ] <- matrix(fit.mat[iter, d4.ind], nr = nt)
    }
  }

  d.mat
}

### Function to deal with under-reporting
### if present
urep.fct <- function(urep, AIDS, urep.flag, qt.flag, nt) { ### including under-reporting in expected AIDS, for each iter
  out <- matrix(NA, nr = nrow(AIDS), nc = ncol(AIDS))
  ### Under-reporting starts from year 2000
  ### which is year 5 (qt.flag=FALSE)
  ### or quarter 20 (qt.flag=FALSE)
  urep.start <- ifelse(qt.flag == TRUE, 20, 5)

  if (urep.flag) {
    out[1:urep.start, ] <- AIDS[1:urep.start, ] ### non under-reporting years
    out[(urep.start + 1):nt, ] <- urep * AIDS[(urep.start + 1):nt, ] ### under-reporting years
  } else {
    out <- AIDS
  }
  out
}
