# library(Rcpp)
# library(Rsrat)
# library(lhs)
#
# sourceCpp(file="lhc.cpp")
# sourceCpp(file="KS.cpp")
# source("kstest.R")
#
# exam <- function(data, models) {
#   result <- fit.srm.nhpp(fault=data, srm.names=models, selection = NULL)
#   p <- mvfplot(time=result$exp$srm$data$time, fault=result$exp$srm$data$fault,
#                mvf=lapply(result, function(m) m$srm))
#
#   ksresult0 <- lapply(result, function(m) {
#     tres <- system.time(res <- ks.srm.test.y(m$srm))
#     c(res, list(comp.time=tres[1]))
#   })
#   names(ksresult0) <- models
#   labels <- names(ksresult0$exp)
#   result.table0 <- lapply(labels, function(l) sapply(ksresult0, function(x) x[[l]]))
#   names(result.table0) <- labels
#   df0 <- data.frame(result.table0)
#   df0 <- transform(df0, p95=df0$p.value < 0.05, p99=df0$p.value<0.01)
#   names(df0) <- sapply(names(df0), function(x) paste("yamada.", x, sep = ""))
#
#   ksresult <- lapply(result, function(m) {
#     tres <- system.time(res <- gks.test(m$srm))
#     c(res, list(comp.time=tres[1]))
#   })
#   names(ksresult) <- models
#   labels <- names(ksresult$exp)
#   result.table <- lapply(labels, function(l) sapply(ksresult, function(x) x[[l]]))
#   names(result.table) <- labels
#   df1 <- data.frame(result.table)
#   pfunc <- function(l, u, alpha) {
#     apply(cbind(l,u), 1, function(x) {
#       if (x[2] < alpha) {
#         TRUE
#       } else if (x[1] > alpha) {
#         FALSE
#       } else {
#         NA
#       }
#     })
#   }
#   df1 <- transform(df1,
#                    p95=pfunc(df1$p.value.lower, df1$p.value.upper, 0.05),
#                    p99=pfunc(df1$p.value.lower, df1$p.value.upper, 0.01))
#   names(df1) <- sapply(names(df1), function(x) paste("gks.", x, sep = ""))
#
#   ksresult.mc <- lapply(result, function(m) {
#     tres <- system.time(res <- gks.test.mc(m$srm))
#     c(res, list(comp.time=tres[1]))
#   })
#   names(ksresult.mc) <- models
#   labels <- names(ksresult.mc$exp)
#   result.table.mc <- lapply(labels, function(l) sapply(ksresult.mc, function(x) x[[l]]))
#   names(result.table.mc) <- labels
#   df2 <- data.frame(result.table.mc)
#   df2 <- transform(df2, p95=df2$p.value < 0.05, p99=df2$p.value<0.01)
#   names(df2) <- sapply(names(df2), function(x) paste("gks.mc.", x, sep = ""))
#
#   ksresult.mc2 <- lapply(result, function(m) {
#     tres <- system.time(res <- gks.test.mc2(m$srm))
#     c(res, list(comp.time=tres[1]))
#   })
#   names(ksresult.mc2) <- models
#   labels <- names(ksresult.mc2$exp)
#   result.table.mc2 <- lapply(labels, function(l) sapply(ksresult.mc2, function(x) x[[l]]))
#   names(result.table.mc2) <- labels
#   df3 <- data.frame(result.table.mc2)
#   df3 <- transform(df3, p95=df3$p.value < 0.05, p99=df3$p.value<0.01)
#   names(df3) <- sapply(names(df3), function(x) paste("gks.mc2.", x, sep = ""))
#
#   list(df=cbind(df0, df1, df2, df3), graph=p)
# }
#
# data(dacs)
# # result.sys17g <- exam(sys17g, srm.models)
# # print(result.sys17g)
# # write.csv(result.sys17g$df, file="sys17g.csv")
# # result.ss3g.168 <- exam(ss3g[1:168], srm.models)
# # print(result.ss3g.168)
# # write.csv(result.ss3g.168$df, file="ss3g168.csv")
# # result.sys1g <- exam(sys1g, srm.models)
# # print(result.sys1g)
# # write.csv(result.sys1g$df, file="sys1g.csv")
# # result.tohma <- exam(tohma, srm.models)
# # print(result.tohma)
# # write.csv(result.tohma$df, file="tohma.csv")
#
# gdata <- list(
#   ss1ag=ss1ag,
#   ss1bg=ss1bg,
#   ss1cg=ss1cg,
#   ss2g=ss2g,
#   ss3g=ss3g,
#   ss4g=ss4g,
#   sys14cg=sys14cg,
#   sys17g=sys17g,
#   sys1g=sys1g,
#   sys27g=sys27g,
#   sys2g=sys2g,
#   sys3g=sys3g,
#   sys40g=sys40g,
#   sys4g=sys4g,
#   sys5g=sys5g,
#   sys6g=sys6g,
#   tohma=tohma
# )
#
# result.all <- lapply(1:length(gdata), function(i) {
#   result <- exam(gdata[[i]], srm.models)
#   write.csv(result$df, file=paste(names(gdata)[[i]], ".csv", sep = ""))
#   result
# })
