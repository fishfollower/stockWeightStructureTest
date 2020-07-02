library(TMB)
compile("lm.cpp")
dyn.load(dynlib("lm"))

data <- list()
data$x <- runif(1000, 0,10)
data$Y <- 1*data$x+2+rnorm(1000)+5

param <- list()
param$alpha <- 0
param$beta <- 0
param$logSigma <- 0
param$logLambda <- 0

obj <- MakeADFun(data, param, DLL="lm")
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))
