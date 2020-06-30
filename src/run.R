library(TMB)
# compile 
compile("../../src/gmrf1.cpp")
dyn.load(dynlib("../../src/gmrf1"))

cat("###################\n",getwd(),"\n#####################")

runit <- function(mode=1, transCode=-1, res=FALSE, label=NULL, predict=0, cut=predict-1, silent=TRUE, cut.data=0, map=list(), ...){
  
  ## trans codes: -1 = identity, 0 = log, >0 = power [ x^(transCode) ]
  trans <- identity
  invtrans <- identity
  lab <- paste0(mode," , identity")  
  if( transCode == 0 ){
      trans <- log
      invtrans <- exp
      lab<-paste0(mode," , log")
  } else if( transCode > 0){
      trans <- function(x) x^(transCode)
      invtrans <- function(x) x^(1/transCode)
      lab <- paste0(mode," , x^",transCode)
  }
  if(missing(label))label<-lab #paste0(mode,",",deparse(substitute(trans)))
    
  # setup data 
  data <- list()

  Yorg <- as.matrix(read.table("Y.tab", head=FALSE))
  Norg <- as.matrix(read.table("N.tab", head=FALSE))
  Moorg <- as.matrix(read.table("Mo.tab", head=FALSE))
  Yorg[abs(Yorg)<1.0e-12] <- NA
  Y <- trans(Yorg)
  if(cut.data>0){
    Yorg<-Yorg[1:(nrow(Yorg)-cut.data),]
    Y<-Y[1:(nrow(Y)-cut.data),]
    Norg<-Norg[1:(nrow(Norg)-cut.data),]
    Moorg<-Moorg[1:(nrow(Moorg)-cut.data),]      
  }    
  Y[!(1:nrow(Y)%in%(1:(nrow(Y)-predict))),]<-NA
  if(cut>0){
    Yorg<-Yorg[1:(nrow(Yorg)-cut),]
    Y<-Y[1:(nrow(Y)-cut),]
    Norg<-Norg[1:(nrow(Norg)-cut),]
    Moorg<-Moorg[1:(nrow(Moorg)-cut),]      
  }
  jac <- -sum(log(abs(numDeriv:::grad(trans,Yorg[!is.na(Y)]))))
  
  r <- as.vector(row(Y))
  c <- as.vector(col(Y))
  n <- length(r)
  W.r <- W.c <- W.d <- matrix(0,nrow=n, ncol=n)

  for(i in 1:n){
    for(j in 1:n){
      W.r[i,j] <- (r[i]==r[j])&(abs(c[i]-c[j])==1)
      W.c[i,j] <- (c[i]==c[j])&(abs(r[i]-r[j])==1)
      W.d[i,j] <- (((r[i]-r[j])==1)&((c[i]-c[j])==1))|(((r[i]-r[j])==(-1))&((c[i]-c[j]) ==(-1)))
    }      
  }
  diag(W.r) <- -rowSums(W.r)
  diag(W.c) <- -rowSums(W.c)
  diag(W.d) <- -rowSums(W.d)

  data$mode <- mode
  data$Wr <- W.r
  data$Wc <- W.c
  data$Wd <- W.d
  data$Y <- Y

  # setup parameters
  param<-list()  
  if(mode==0){
    data$aveYears <- 3      
    param$logSdObs <- numeric(ncol(data$Y))
    ran <- NULL
  }
  if(mode==1){    
    param$logPhi <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    ran <- c("z")
  }
  if(mode==2){    
    param$logitRho <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- c(0,0)
    param$logSdObs <- numeric(ncol(data$Y))
    param$omega <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$z <- rep(0,(nrow(Y)-1)+ncol(Y))
    ran <- c("omega","z")
  }
  if(mode==3 || mode==4){
    param$logitRho <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- c(0,0)
    param$logSdObs <- numeric(ncol(data$Y))
    param$omega <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$z <- rep(0,(nrow(Y)-1)+ncol(Y))
    param$logitRhoObs <- 0
    if(mode==4){
        data$trans <- transCode
        stopifnot(transCode>=0)
    }
    ran <- c("omega","z")
  }
  

  ## run model 
  obj <- MakeADFun(data,param,random=ran, DLL="gmrf1", silent=silent, map=map, ...)
    
  lower <- rep(-Inf,length(obj$par))
  upper <- rep(Inf,length(obj$par))  
  pn <- names(obj$par)
  lower[ grep("^log",pn) ] <- -5
  upper[ grep("^log",pn) ] <- 5
  lower[ grep("^logit",pn) ] <- -4
  upper[ grep("^logit",pn) ] <- 4
    
    
  opt <- nlminb(obj$par, obj$fn, obj$gr,lower=lower,upper=upper)
  sdr <- sdreport(obj)
  pred <- as.list(sdr, report=TRUE, what="Est")$pred
  predSd <- as.list(sdr, report=TRUE, what="Std")$pred

  sink("sdrep.tab",append=TRUE); print(label); print(summary.sdreport(sdr,"fixed")); sink();
  loglik <- opt$objective+jac
  k <- length(opt$par)
  nobs <- length(data$Y)  
  aicc <- 2*loglik + 2*k + 2*k*(k+1)/(nobs-k-1)
  aic <- 2*loglik + 2*k
    
  matplot(invtrans(pred), type="l", ylim=range(invtrans(data$Y), na.rm=TRUE), main=label)
  matplot(invtrans(pred-2*predSd), , type="l", add=TRUE, lty="dotted")
  matplot(invtrans(pred+2*predSd), type="l", add=TRUE, lty="dotted")
  matplot(invtrans(data$Y), add=TRUE)
    
  residual <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))     
  if(res){
    ooa <- oneStepPredict(obj, data.term.indicator="keep", observation.name="Y", discrete=FALSE)
    residual <- matrix(ooa$residual, nrow=nrow(Y), ncol=ncol(Y))
    stockassessment:::plotby(row(residual), col(residual), residual)
    boxplot(residual)   
  }
    
  SSBorg<-rowSums(Yorg*Norg*Moorg)
  SSBpred<-rowSums(invtrans(pred)*Norg*Moorg)    
  conv<-(all(is.finite(summary.sdreport(sdr, "fixed")[,2])))&(opt$convergence==0)
  if(!conv)warning("Convergence issue")  
  return(list(logLik=opt$objective+jac, AICc=aicc, AIC=aic, obj=obj, residual=residual, opt=opt, sdr=sdr, ssbobs=SSBorg, ssbpred=SSBpred, conv=conv, label=label, call=mget(names(formals()),sys.frame(sys.nframe()))))
}

cv.rmse <- function(year=10, cv.scale=identity, lag=10, ...){
    sq.error<-function(p){fit<-runit(...,predict=p); c((tail(cv.scale(fit$ssbobs),1)-tail(cv.scale(fit$ssbpred),1))^2,as.integer(fit$conv))}
    ret<-rowMeans(Vectorize(sq.error)(1:year+lag))
    ret[1]<-sqrt(ret[1])
    ret
}

dat <- read.table("Y.tab", head=FALSE)
mymap <- list(logSdObs=factor(rep(1,ncol(dat))))

pdf("res.pdf")
  mod <- list()
  mod[[length(mod)+1]] <- runit(mode=0, trans=0, res=TRUE, map=mymap, cut.data=10, label="Mod0-log-constantVariance")
  ##mod[[length(mod)+1]] <- runit(mode=1, res=TRUE, map=mymap, cut.data=10, label="Mod1-identity-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=1, trans=0, res=TRUE, map=mymap, cut.data=10, label="Mod1-log-constantVariance")
  ##mod[[length(mod)+1]] <- runit(mode=2, res=TRUE, map=mymap, cut.data=10, label="Mod2-identity-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=2, trans=0, res=TRUE, map=mymap, cut.data=10, label="Mod2-log-constantVariance")
  ##mod[[length(mod)+1]] <- runit(mode=3, res=TRUE, map=mymap, cut.data=10, label="Mod3-identity-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=3, trans=0, res=TRUE, map=mymap, cut.data=10, label="Mod3-log-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=4, trans=0, res=TRUE,map=mymap, cut.data=10, label="Mod4-log-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=4, trans=1/3, res=TRUE,map=mymap, cut.data=10, label="Mod4-cubrt-constantVariance")
  mod[[length(mod)+1]] <- runit(mode=4, trans=1/2, res=TRUE,map=mymap, cut.data=10, label="Mod4-sqrt-constantVariance")
dev.off()

res <- as.data.frame(do.call(rbind, lapply(mod, function(m)c(m$label, round(m$logLik,2), round(m$AICc,2), round(m$AIC,2), m$conv))))

cv<-lapply(mod, function(m)cv.rmse(year=10, cv.scale=log, mode=m$call$mode, transCode=m$call$transCode, label=m$label, cut.data=10, map=m$call$map))

res<-cbind(res,do.call(rbind,cv))

names(res)<-c("Label", "-logLik", "AICc","AIC","Conv all", "RMSE-CV", "Conv rate CV")

cat(sub("^[1-9]*","",capture.output(res)), file = 'res.tab', sep = '\n')

