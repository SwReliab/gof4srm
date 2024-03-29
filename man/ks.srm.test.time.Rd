% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kstest.R
\name{ks.srm.test.time}
\alias{ks.srm.test.time}
\title{KS test for software reliability models with time data}
\usage{
ks.srm.test.time(obj, alternative = c("two.sided", "less", "greater"))
}
\arguments{
\item{obj}{An object of an estimated result by Rsrat}

\item{alternative}{A string indicates the alternative hypothesis and
must be one of "two.sided" (default), "less", or "greater".}
}
\value{
A list with components;
\item{statistic}{A value of the test statistic.}
\item{p.value}{A value of the p-value of the test.}
\item{alternative}{A string of the alternative hypothesis.}
}
\description{
Perform the KS test for the estimated software reliability models
}
\examples{
data(dacs)
result <- fit.srm.nhpp(time=sys1[sys1>=0],
                       te=-sys1[sys1<0], srm.name=c("exp"))
ks.srm.test.time(result)
}
