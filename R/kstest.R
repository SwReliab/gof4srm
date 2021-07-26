#' Truncated CDF for software reliability model
#'
#' Generate the truncated CDF from the estiamted software
#' reliability model
#'
#' @param srm An object for the estimated software reliability model
#' @return A function object for the truncated CDF.

trunc.cdf <- function(srm) {
  function(t) {
    te <- sum(srm$data$time)
    srm$mvf(t) / srm$mvf(te)
  }
}

#' Resampling for the grouped data
#'
#' @param n An integer for the number of samples
#' @param size An integer indicates the number of total faults
#' @param ctime A sequence for time slots
#' @param cdf A function object of the truncated CDF
#' @return A matrix of resamples. Each column corresponds to a sample

resample <- function(n, size, ctime, cdf) {
  p <- diff(c(0, cdf(ctime)))
  rmultinom(n, size=size, prob=p)
}

#' Resampling for the grouped data with Latin hypercube
#'
#' @param n An integer for the number of samples
#' @param size An integer indicates the number of total faults
#' @param ctime A sequence for time slots
#' @param cdf A function object of the truncated CDF
#' @return A matrix of resamples. Each column corresponds to a sample

resample.LHC <- function(n, size, ctime, cdf) {
  multinomialLHC(n, size, ctime, cdf, randomLHS(n, size))
}

#' KS test for software reliability models
#'
#' Perform the KS test for the estimated software reliability models
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @param method A string indicates the method for KS test with group data. "mc"
#' and "mc2" are Monte-Carlo methods. "gks" is the generalized KS test. "yamada"
#' is an approximation. The default is "mc".
#' @param ... A list of other options that are passed to concrete functions.
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value}{A value of the p-value of the test.}
#' \item{p.value.lower}{A lower bound of the p-value of the test. (gks)}
#' \item{p.value.upper}{A upper bound of the p-value of the test. (gks)}
#' \item{alternative}{A string of the alternative hypothesis.}
#' \item{b}{An integer for the number of resamples. (mc, mc2)}
#' \item{alpha}{A value of the significant level. (mc, mc2)}
#' \item{cv}{A value of the coefficient variation. (mc, mc2)}
#' \item{lhc}{A logical indicates resamples are drawn by Latin hypercude.
#' (mc, mc2)}
#' \item{seed}{An integer for the seed of random numbers. (mc, mc2)}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(time=sys1[sys1>=0],
#'                        te=-sys1[sys1<0], srm.name=c("exp"))
#' ks.srm.test(result)
#' results <- fit.srm.nhpp(time=sys1[sys1>=0],
#'                         te=-sys1[sys1<0], srm.name=c("exp", "gamma"))
#' sapply(results, function(x) x$p.value)
#' @export

ks.srm.test <- function(obj, alternative = c("two.sided", "less", "greater"),
                        method = c("mc", "gks", "mc2", "yamada"), ...) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  if (check.faultdata(obj$srm$data) == "time") {
    result <- ks.srm.test.time(obj, alternative = alternative)
    c(result, list(method="ks"))
  } else if (check.faultdata(obj$srm$data) == "group") {
    result <- switch(method,
      "gks"    = gks.srm.test.mp(obj, alternative = alternative, ...),
      "mc" = gks.srm.test.mc(obj, alternative =
                               c("two.sided", "less", "greater"), ...),
      "mc2" = gks.srm.test.mc2(obj, alternative =
                                 c("two.sided", "less", "greater"), ...),
      "yamada" = ks.srm.test.y(obj, alternative = alternative),
      stop("Method is not defined")
    )
    c(result, list(method=method))
  } else {
    stop("The fault data is neither time nor group data.")
  }
}

#' KS test for software reliability models with time data
#'
#' Perform the KS test for the estimated software reliability models
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value}{A value of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(time=sys1[sys1>=0],
#'                        te=-sys1[sys1<0], srm.name=c("exp"))
#' ks.srm.test.time(result)
#' @export

ks.srm.test.time <- function(obj, alternative =
                               c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  if (alternative != "two.sided") {
    stop("less and greater are not implemented yet")
  }
  tcdf <- trunc.cdf(obj$srm)
  ctime <- cumsum(obj$srm$data$time)
  fault <- obj$srm$data$type
  size <- sum(fault)

  ks0 <- KSdistance_point(sample=ctime, tcdf)
  statistic <- switch(alternative,
                      "two.sided" = ks0$d,
                      "less" =  ks0$dminus,
                      "greater" = ks0$dplus)
  pvalue <- 1 - .Call(stats:::C_pKolmogorov2x, STAT=statistic, n=size)
  list(statistic=statistic, p.value=pvalue, alternative=alternative)
}

#' KS test for software reliability models with group data
#' (Yamada's approximation)
#'
#' Perform the KS test for the estimated software reliability models
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value}{A value of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=sys1g[1:30], srm.name=c("exp"))
#' ks.srm.test.y(result)
#' @export

ks.srm.test.y <- function(obj, alternative =
                            c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  if (alternative != "two.sided") {
    stop("less and greater are not implemented yet")
  }
  tcdf <- trunc.cdf(obj$srm)
  ctime <- cumsum(obj$srm$data$time)
  fault <- obj$srm$data$fault
  size <- sum(fault)

  ks0 <- KSdistance_group2(ctime=ctime, count=fault, tcdf)
  statistic <- switch(alternative,
                      "two.sided" = ks0$d,
                      "less" =  ks0$dminus,
                      "greater" = ks0$dplus)
  pvalue <- 1 - .Call(stats:::C_pKolmogorov2x, STAT=statistic, n=size)
  list(statistic=statistic, p.value=pvalue, alternative=alternative)
}

#' Generalized KS test for software reliability models with group data
#'
#' Perform the generalized KS test for the estimated software reliability models
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value.lower}{A lower bound of the p-value of the test.}
#' \item{p.value.upper}{A upper bound of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}

gks.srm.test <- function(obj, alternative = c("two.sided", "less", "greater")) {
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
  prob.dm <- compute_Pdminus(ctime=ctime, count=fault, statistic, tcdf)
  prob.dp <-compute_Pdplus(ctime=ctime, count=fault, statistic, tcdf)

  pvalue.lower <- prob.dp + prob.dm - prob.dp * prob.dm
  pvalue.upper <- prob.dp + prob.dm
  list(statistic=statistic, p.value.lower=pvalue.lower,
       p.value.upper=pvalue.upper, alternative=alternative)
}

#' Generalized KS test for software reliability models with Monte-Carlo testing
#'
#' Perform the generalized KS test for the estimated software reliability models
#' with Monte-Carlo testing
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @param b An integer for the number of resamples. If it is NULL, the number of resamples
#' is determined from the significant level and the coefficient variation.
#' @param alpha A value of the significant level which is used to determine the number
#' of resamples. The default is 0.01. If \code{b} is not NULL, this is ignored.
#' @param cv A value of the coefficient variation which is used to determine the number
#' of resamples. The default is 0.1. If \code{b} is not NULL, this is ignored.
#' @param lhc A logical indicates resamples are drawn by Latin hypercude.
#' @param seed An integer for the seed of random numbers. If it is NULL, the seed is not set.
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value}{A value of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}
#' \item{b}{An integer for the number of resamples.}
#' \item{alpha}{A value of the significant level.}
#' \item{cv}{A value of the coefficient variation.}
#' \item{lhc}{A logical indicates resamples are drawn by Latin hypercude.}
#' \item{seed}{An integer for the seed of random numbers.}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=sys1g[1:30], srm.name=c("exp"))
#' gks.srm.test.mc(result)
#' @export

gks.srm.test.mc <- function(obj, alternative = c("two.sided", "less", "greater"),
                      b = NA, alpha = 0.01, cv = 0.1, lhc = FALSE, seed = NA) {
  alternative <- match.arg(alternative)
  if (is.na(b)) {
    b <- ceiling((1-alpha) / (cv^2 * alpha))
  }
  tcdf <- trunc.cdf(obj$srm)
  ctime <- cumsum(obj$srm$data$time)
  fault <- obj$srm$data$fault
  size <- sum(fault)

  ks0 <- KSdistance_group(ctime=ctime, count=fault, tcdf)
  statistic <- switch(alternative,
                      "two.sided" = ks0$d,
                      "less" =  ks0$dminus,
                      "greater" = ks0$dplus)

  if (!is.na(seed)) {
    set.seed(seed)
  }
  if (lhc) {
    sample <- resample.LHC(b, size, ctime, tcdf)
  } else {
    sample <- resample(b, size, ctime, tcdf)
  }
  res <- KSdistance_groupM(ctime, size, sample, tcdf);
  p.value <- switch(alternative,
                    "two.sided" = (1+sum(res$d >= statistic)) / (b+1),
                    "less" =  (1+sum(res$dminus >= statistic)) / (b+1),
                    "greater" = (1+sum(res$dplus >= statistic)) / (b+1))

  list(statistic=statistic, p.value=p.value, alternative=alternative,
       b=b, alpha=alpha, cv=cv, lhc=lhc, seed=seed)
}

#' Generalized KS test for software reliability models with Monte-Carlo testing
#' under unknown model parameters
#'
#' Perform the generalized KS test for the estimated software reliability models
#' with Monte-Carlo testing
#'
#' @param obj An object of an estimated result by Rsrat
#' @param alternative A string indicates the alternative hypothesis and
#' must be one of "two.sided" (default), "less", or "greater".
#' @param b An integer for the number of resamples. If it is NULL, the number of resamples
#' is determined from the significant level and the coefficient variation.
#' @param alpha A value of the significant level which is used to determine the number
#' of resamples. The default is 0.01. If \code{b} is not NULL, this is ignored.
#' @param cv A value of the coefficient variation which is used to determine the number
#' of resamples. The default is 0.1. If \code{b} is not NULL, this is ignored.
#' @param lhc A logical indicates resamples are drawn by Latin hypercude.
#' @param seed An integer for the seed of random numbers. If it is NULL, the seed is not set.
#' @return A list with components;
#' \item{statistic}{A value of the test statistic.}
#' \item{p.value}{A value of the p-value of the test.}
#' \item{alternative}{A string of the alternative hypothesis.}
#' \item{b}{An integer for the number of resamples.}
#' \item{alpha}{A value of the significant level.}
#' \item{cv}{A value of the coefficient variation.}
#' \item{lhc}{A logical indicates resamples are drawn by Latin hypercude.}
#' \item{seed}{An integer for the seed of random numbers.}
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=sys1g[1:30], srm.name=c("exp"))
#' gks.srm.test.mc2(result)
#' @export

gks.srm.test.mc2 <- function(obj, alternative = c("two.sided", "less", "greater"),
                        b = NA, alpha = 0.01, cv = 0.1, lhc = FALSE, seed = NA) {
  alternative <- match.arg(alternative)
  if (is.na(b)) {
    b <- ceiling((1-alpha) / (cv^2 * alpha))
  }
  tcdf <- trunc.cdf(obj$srm)
  ctime <- cumsum(obj$srm$data$time)
  fault <- obj$srm$data$fault
  size <- sum(fault)

  ks0 <- KSdistance_group(ctime=ctime, count=fault, tcdf)
  statistic <- switch(alternative,
                      "two.sided" = ks0$d,
                      "less" =  ks0$dminus,
                      "greater" = ks0$dplus)

  if (!is.na(seed)) {
    set.seed(seed)
  }
  if (lhc) {
    ss <- resample.LHC(b, size, ctime, tcdf)
  } else {
    ss <- resample(b, size, ctime, tcdf)
  }
  data <- obj$srm$data
  dvec <- apply(ss, 2, function(x) {
    data$fault <- x
    ressrm <- emfit(srm, data)$srm
    res <- KSdistance_group(ctime, x, trunc.cdf(ressrm))
    switch(alternative,
           "two.sided" = res$d,
           "less" =  res$dminus,
           "greater" = res$dplus)
  })
  p.value <- (1+sum(dvec >= statistic)) / (b+1)
  list(statistic=statistic, p.value=p.value, alternative=alternative,
       b=b, alpha=alpha, cv=cv, lhc=lhc, seed=seed)
}
