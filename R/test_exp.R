# # test file for the ordinary KS for exponential distirbution
# # We compared the p-values
#
# lambda <- 0.1
# n <- 50
# x <- rexp(n, lambda)
# #x <- rgamma(n, shape=1.5, rate=1)
#
# n <- length(x)
# s <- sum(x)
#
# dn <- function(x, lambda) {
#   n <- length(x)
#   max(apply(
#   cbind(
#   sapply(sort(x), function(u) pexp(u, lambda)),
#   (1:n - 1)/n,
#   1:n / n), 1, function(v) max(abs(v[1]-v[2]), abs(v[1]-v[3]))))
# }
# ks0 <- dn(x, n/s)
#
# bs <- c()
# bs1 <- c()
# bs2 <- c()
# for (i in 1:10000) {
#   b <- diff(c(0, sort(runif(n = n-1, min=0, max=s)), s)) #pattern1
#   bs <- c(bs, dn(b, n/s))
#   b1 <- rexp(n, n/s) #pattern2
#   bs1 <- c(bs1, dn(b1, n/s))
#   bs2 <- c(bs2, dn(b1, n/sum(b1)))
# }
#
#
# q0<-quantile(bs,c(1-0.2,1-0.1,1-0.05,1-0.02,1-0.01))
# q1<-quantile(bs1,c(1-0.2,1-0.1,1-0.05,1-0.02,1-0.01))
# q2<-quantile(bs2,c(1-0.2,1-0.1,1-0.05,1-0.02,1-0.01))
#
#
# ks.test(x=x,y=pexp,rate=n/s)
#
# hist(bs)
# pvalue0 <- (1+sum(bs>=ks0))/(10000+1) #MonteCarlo p value
# hist(bs1)
# pvalue1 <- (1+sum(bs1>=ks0))/(10000+1) #MonteCarlo p value
# hist(bs2)
# pvalue2 <- (1+sum(bs2>=ks0))/(10000+1) #MonteCarlo p value
#
#
# pvalue0
# pvalue2
