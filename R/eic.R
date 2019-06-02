# testlam <- function(x, ...) {
#   args <- list(...)
#   r <- ifelse(is.null(args$rate), 1, args$rate)
#   return(args$total * dexp(x, rate=r))
# }
#
# testlam2 <- function(x, ...) {
#   args <- list(...)
#   l <- ifelse(is.null(args$location), 0, args$location)
#   s <- ifelse(is.null(args$scale), 0, args$scale)
#   return(args$total * dlogis(x, location=l, scale=s))
# }

#' Generate samples drawn from NHPP
#'
#' This is a function to generate samples drwan from NHPP by thining.
#'
#' @param intensity A function of intensity function
#' @param range A vector (min, max) indicating the range
#' @param maximum A value of maximum intensity.
#' @return A time series of event occurences

thinning <- function(intensity, range, maximum) {
  n <- rpois(1, maximum*(range[2] - range[1]))
  candidates <- sort(runif(n, min=range[1], max=range[2]))
  y <- candidates[intensity(candidates) / maximum > runif(n)]
  return(sort(y))
}

#' bootstrap of NHPP
#'
#' Generate time data (bootstrap samples) from a given NHPP.
#'
#' @param srm An object of an estimated result by Rsrat
#' @return A time data (time series)
#' @export

bs.nhpp.time <- function(intensity, range, maximum) {
  x <- thinning(intensity, range, maximum)
  time <- diff(c(0, x))
  faultdata(time=time, te=range[2] - sum(time))
}

#' EIC for SRM
#'
#' Compute EIC for SRM
#'
#' @param srm An object of an estimated result by Rsrat
#' @return resA time data (time series)
#' @export

eic <- function(result, bsample = 100, alpha = 0.95) {
  b3 <- result$llf
  b1 <- numeric(bsample)
  b2 <- numeric(bsample)
  b4 <- numeric(bsample)

  range <- c(0, result$srm$data$max)
  maximum <- optimize(result$srm$intensity, range, maximum=TRUE)$objective

  for (b in 1:bsample) {
    sample <- bs.nhpp.time(result$srm$intensity, range, maximum)
    b2[b] <- result$srm$llf(sample)
    result.bs <- emfit(result$srm$clone(), sample, initialize=FALSE)
    b1[b] <- result.bs$llf
    b4[b] <- result.bs$srm$llf(result$srm$data)
  }
  x <- mean(b1 - b2 + b3 - b4)
  s <- stats::sd(b1 - b2 + b3 - b4)
  talpha <- stats::qt(p=(1-alpha)/2, df=bsample-1, lower.tail=FALSE)
  x.interval <- x + c(-talpha, talpha)*s/sqrt(bsample)
  list(bias = x,
       bias.lower = x.interval[1],
       bias.upper = x.interval[2],
       eic = -2 * (result$llf - x),
       eic.lower = -2 * (result$llf -  x.interval[1]),
       eic.upper = -2 * (result$llf -  x.interval[2]))
}
