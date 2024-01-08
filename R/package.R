#' gof4srm: Goodness-of-fit tests for software reliability models
#'
#' This is a package of goodness-of-fit tests for software reliability (growth)
#' models. This provides Kolmogorov-Smirnov test for both time and grouped data,
#' and EIC (extended information criterion).This package provides estimation
#' programs for software reliability growth
#'
#' @docType package
#' @name gof4srm
#' @import Rsrat
#' @import foreach
#' @import parallel
#' @import doParallel
#' @importFrom Rmpfr mpfr
#' @importFrom stats rpois runif rmultinom
#' @importFrom lhs randomLHS
#' @useDynLib gof4srm
NULL
