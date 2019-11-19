#' bootstrap of NHPP
#'
#' Generate time data (bootstrap samples) from a given NHPP.
#'
#' @param srm An object of model provided from Rsrat
#' @param data An object of faultdata
#' @param maximum A value of the maximum value of intensity for an interval (0, the last observation time)
#' @return An object of faultdata of bootstrap sample
#' @export

bs.srm.time <- function(srm, data, maximum) {
  te <- data$max
  cte <- sum(data$time) + te
  n <- rpois(1, maximum*cte)
  candidates <- sort(runif(n, min=0, max=cte))
  y <- candidates[srm$intensity(candidates) / maximum > runif(n)]
  time <- diff(c(0, y))
  faultdata(time=time, te=cte-sum(time))
}

#' bootstrap of NHPP
#'
#' Generate group data (bootstrap samples) from a given NHPP.
#'
#' @param srm An object of model provided from Rsrat
#' @param data An object of faultdata
#' @return An object of faultdata of bootstrap sample
#' @export

bs.srm.group <- function(srm, data) {
  ctime <- c(0, cumsum(srm$data$time))
  lambda <- diff(srm$mvf(ctime))
  faultdata(fault=sapply(lambda, function(r) rpois(1, lambda = r)))
}

#' EIC for SRM
#'
#' Generate samples of bias for SRM with group data from bootstrap sampling
#'
#' @param obj An object of an estimated result by Rsrat
#' @param bsample An integer of the number of bootstrap samples
#' @param initialize A logical. If it is TRUE, the parameters are reset by initial guesses.
#' The default is FALSE.
#' @return A numeric vector of samples of bias
#' @export

eic.srm.group <- function(obj, bsample = 100, initialize = FALSE) {
  b3 <- obj$llf
  b1 <- numeric(bsample)
  b2 <- numeric(bsample)
  b4 <- numeric(bsample)

  ctime <- c(0, cumsum(obj$srm$data$time))
  lambda <- diff(obj$srm$mvf(ctime))

  for (b in 1:bsample) {
    sample <- faultdata(fault=sapply(lambda, function(r) rpois(1, lambda = r)))
    b2[b] <- obj$srm$llf(sample)
    obj.bs <- emfit(obj$srm$clone(), sample, initialize=initialize)
    b1[b] <- obj.bs$llf
    b4[b] <- obj.bs$srm$llf(obj$srm$data)
  }
  b1 - b2 + b3 - b4
}

#' EIC for SRM
#'
#' Generate samples of bias for SRM with time data from bootstrap sampling
#'
#' @param obj An object of an estimated result by Rsrat
#' @param bsample An integer of the number of bootstrap samples
#' @param initialize A logical. If it is TRUE, the parameters are reset by initial guesses.
#' The default is FALSE.
#' @return A numeric vector of samples of bias
#' @export

eic.srm.time <- function(obj, bsample = 100, initialize = FALSE) {
  b3 <- obj$llf
  b1 <- numeric(bsample)
  b2 <- numeric(bsample)
  b4 <- numeric(bsample)

  cte <- sum(obj$data$time) + obj$srm$data$max
  maximum <- stats::optimize(obj$srm$intensity, c(0,cte), maximum=TRUE)$objective

  for (b in 1:bsample) {
    sample <- bs.srm.time(obj$srm, obj$srm$data, maximum)
    b2[b] <- obj$srm$llf(sample)
    obj.bs <- emfit(obj$srm$clone(), sample, initialize=initialize)
    b1[b] <- obj.bs$llf
    b4[b] <- obj.bs$srm$llf(obj$srm$data)
  }
  b1 - b2 + b3 - b4
}

#' EIC for SRM
#'
#' Generate samples of EIC
#'
#' @param obj An object of an estimated result by Rsrat
#' @param bsample An integer of the number of bootstrap samples
#' @param initialize A logical. If it is TRUE, the parameters are reset by initial guesses.
#' The default is FALSE.
#' @return A vector of samples
#' @export

eic.srm.sample <- function(obj, bsample = 100, initialize = FALSE) {
  if (check.faultdata(obj$srm$data) == "time") {
    bias.sample <- eic.srm.time(obj, bsample, initialize)
  } else if (check.faultdata(obj$srm$data) == "group") {
    bias.sample <- eic.srm.group(obj, bsample, initialize)
  } else {
    stop("The fault data is neither time nor group data.")
  }
  -2 * (obj$llf - bias.sample)
}

#' EIC for SRM
#'
#' Compute EIC for SRM
#'
#' @param obj An object of an estimated result by Rsrat
#' @param bsample An integer of the number of bootstrap samples
#' @param alpha A value of significant level to obtain confidence interval.
#' The default is 0.95.
#' @param initialize A logical. If it is TRUE, the parameters are reset by initial guesses.
#' The default is FALSE.
#' @return lower and upper bounds of bias, eic
#' @export

eic.srm <- function(obj, bsample = 100, alpha = 0.95, initialize = FALSE) {
  if (check.faultdata(obj$srm$data) == "time") {
    bias.sample <- eic.srm.time(obj, bsample, initialize)
  } else if (check.faultdata(obj$srm$data) == "group") {
    bias.sample <- eic.srm.group(obj, bsample, initialize)
  } else {
    stop("The fault data is neither time nor group data.")
  }
  x <- mean(bias.sample)
  s <- stats::sd(bias.sample)
  talpha <- stats::qt(p=(1-alpha)/2, df=bsample-1, lower.tail=FALSE)
  x.interval <- x + c(-talpha, talpha)*s/sqrt(bsample)
  list(bias = x,
       bias.lower = x.interval[1],
       bias.upper = x.interval[2],
       eic = -2 * (obj$llf - x),
       eic.lower = -2 * (obj$llf -  x.interval[1]),
       eic.upper = -2 * (obj$llf -  x.interval[2]))
}
