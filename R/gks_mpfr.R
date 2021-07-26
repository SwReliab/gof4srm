#' Compute P(D->d-) for the generalized KS statistic
#'
#' @param ctime A sequence represents time slots (bins)
#' @param count A sequence indicates the number of samples falls int a bin
#' @param dminus A value of d-
#' @param cdf A function of CDF. It is allowd to get a closure.
#' @param prec An integer to indicate the precision
#' @return A value of the probability
#' @export

comp.Pdminus.mp <- function(ctime, count, dminus, cdf, prec=128) {
  K <- length(ctime)
  n <- sum(count)
  jmax <- floor(n * (1 - dminus))

  cv <- Rmpfr::mpfr(numeric(jmax+1),prec)
  cp <- Rmpfr::mpfr(cdf(ctime),prec)
  for (j in 0:jmax) {
    p <- Rmpfr::mpfr(dminus,prec) + Rmpfr::mpfr(j,prec) / Rmpfr::mpfr(n,prec)
    for (i in 0:K) {
      if (cp[i+1] >= p) {
        cv[j+1] <- Rmpfr::mpfr(1,prec) - cp[i+1]
        break
      }
    }
  }

  bv <- Rmpfr::mpfr(numeric(jmax+1), prec)
  bv[1] <- Rmpfr::mpfr(1,prec)
  for (k in 1:jmax) {
    prob <- Rmpfr::mpfr(0,prec)
    binom <- Rmpfr::mpfr(1,prec)
    for (j in 0:(k-1)) {
      prob <- prob + binom * cv[j+1]^(k-j) * bv[j+1]
      binom <- binom * Rmpfr::mpfr(k-j,prec) / Rmpfr::mpfr(j+1,prec)
    }
    bv[k+1] <- Rmpfr::mpfr(1,prec) - prob
  }

  prob <- Rmpfr::mpfr(0,prec)
  binom <- Rmpfr::mpfr(1,prec)
  for (j in 0:jmax) {
    prob <- prob + binom * cv[j+1]^(n-j) * bv[j+1]
    binom <- binom * Rmpfr::mpfr(n-j,prec) / Rmpfr::mpfr(j+1,prec)
  }

  as.numeric(prob)
}

#' Compute P(D+>d+) for the generalized KS statistic
#'
#' @param ctime A sequence represents time slots (bins)
#' @param count A sequence indicates the number of samples falls int a bin
#' @param dplus A value of d+
#' @param cdf A function of CDF. It is allowd to get a closure.
#' @param prec An integer to indicate the precision
#' @return A value of the probability
#' @export

comp.Pdplus.mp <- function(ctime, count, dplus, cdf, prec=128) {
  K <- length(ctime)
  n <- sum(count)
  jmax <- floor(n * (1 - dplus))

  cv <- Rmpfr::mpfr(numeric(jmax+1),prec)
  cp <- Rmpfr::mpfr(cdf(ctime),prec)
  cp0 <- Rmpfr::mpfr(0,prec)
  for (j in 0:jmax) {
    p <- Rmpfr::mpfr(1,prec) - Rmpfr::mpfr(dplus,prec) - Rmpfr::mpfr(j,prec) / Rmpfr::mpfr(n,prec)
    for (i in 0:K) {
      if (cp[i+1] >= p) {
        cv[j+1] <- cp0
        break
      }
      cp0 <- cp[i+1]
    }
  }

  bv <- Rmpfr::mpfr(numeric(jmax+1),prec)
  bv[1] <- Rmpfr::mpfr(1,prec)
  for (k in 1:jmax) {
    prob <- Rmpfr::mpfr(0,prec)
    binom <- Rmpfr::mpfr(1,prec)
    for (j in 0:(k-1)) {
      prob <- prob + binom * cv[j+1]^(k-j) * bv[j+1]
      binom <- binom * Rmpfr::mpfr(k-j,prec) / Rmpfr::mpfr(j+1,prec)
    }
    bv[k+1] <- Rmpfr::mpfr(1,prec) - prob
  }

  prob <- Rmpfr::mpfr(0,prec)
  binom <- Rmpfr::mpfr(1,prec)
  for (j in 0:jmax) {
    prob <- prob + binom * cv[j+1]^(n-j) * bv[j+1]
    binom <- binom * Rmpfr::mpfr(n-j,prec) / Rmpfr::mpfr(j+1,prec)
  }

  as.numeric(prob)
}

#' Generalized KS test for software reliability models
#'
#' Peform the generalized KS test for the estimated software reliability models
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @param prec An integer to indicate the precision
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value.lower}{A lower bound of the p-value of the test.}
#' \item{p.value.upper}{A upper bound of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=sys1g[1:30], srm.name=c("exp"))
#' gks.srm.test.mp(result)
#' @export

gks.srm.test.mp <- function(obj, alternative = c("two.sided", "less", "greater"), prec = 128) {
  alternative <- match.arg(alternative)
  tcdf <- trunc.cdf(obj$srm)
  ctime <- cumsum(obj$srm$data$time)
  fault <- obj$srm$data$fault
  size <- sum(fault)

  ks0 <- KSdistance_group(ctime=ctime, count=fault, tcdf)
  statistic <- switch(alternative,
                      "two.sided" = ks0$d,
                      "less" =  ks0$dminus,
                      "greater" = ks0$dplus)
  prob.dm <- comp.Pdminus.mp(ctime=ctime, count=fault, statistic, tcdf, prec)
  prob.dp <-comp.Pdplus.mp(ctime=ctime, count=fault, statistic, tcdf, prec)

  pvalue.lower <- prob.dp + prob.dm - prob.dp * prob.dm
  pvalue.upper <- prob.dp + prob.dm
  list(statistic=statistic, p.value.lower=pvalue.lower,
       p.value.upper=pvalue.upper, alternative=alternative)
}
