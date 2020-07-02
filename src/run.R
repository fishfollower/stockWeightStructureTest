library(TMB)
# compile 
compile("../../src/gmrf1.cpp")
dyn.load(dynlib("../../src/gmrf1"))

cat("###################\n",getwd(),"\n#####################")

runit <- function(mode=1, transCode=-1, res=FALSE, label=NULL, predict=0, cut=predict-1, silent=TRUE, cut.data=0, map=list(), lowerLog=-Inf, upperLog=Inf, lowerLogit=-Inf, upperLogit=Inf,  ...){
  
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
  if(mode%in%c(1001,1010,1100,1110,1101,1011,1111)){
    data$mode <- 1
    idx<-as.numeric(strsplit(as.character(mode), "")[[1]])[-1]
    param$logPhi <- c(0,0,0)
    param$logPhi[idx==0]<- -10
    addmap<-1:3
    addmap[idx==0] <- NA
    map=c(map,list(logPhi=factor(addmap)))
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
  if(mode==5){
    param$logPhi <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$logitRhoObs <- 0
    ran <- c("z")  
  }
  if(mode==6){
    param$logPhi <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$logitRhoObs <- 0
    ran <- c("z")
    data$trans <- transCode
    stopifnot(transCode>=0)
  }

  ## run model 
  obj <- MakeADFun(data,param,random=ran, DLL="gmrf1", silent=silent, map=map, ...)
  lower <- rep(-Inf,length(obj$par))
  upper <- rep(Inf,length(obj$par))  
  pn <- names(obj$par)
  lower[ grep("^log",pn) ] <- lowerLog
  upper[ grep("^log",pn) ] <- upperLog
  lower[ grep("^logit",pn) ] <- lowerLogit
  upper[ grep("^logit",pn) ] <- upperLogit
    
    
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

cv.rmse <- function(year=10, cv.scale=identity, ...){
    sq.error<-function(p){fit<-runit(...,predict=p); c((tail(cv.scale(fit$ssbobs),1)-tail(cv.scale(fit$ssbpred),1))^2,as.integer(fit$conv))}
    ret<-rowMeans(Vectorize(sq.error)(1:year))
    ret[1]<-sqrt(ret[1])
    ret
}

jitfun <- function(fit,n=10){
  lower <- rep(-Inf,length(fit$obj$par))
  upper <- rep(Inf,length(fit$obj$par))  
  pn <- names(fit$obj$par)
  lower[ grep("^log",pn) ] <- fit$call$lowerLog
  upper[ grep("^log",pn) ] <- fit$call$upperLog
  lower[ grep("^logit",pn) ] <- fit$call$lowerLogit
  upper[ grep("^logit",pn) ] <- fit$call$upperLogit
  fits <- lapply(1:n, function(i)nlminb(fit$obj$par+rnorm(length(fit$obj$par),sd=.25), fit$obj$fn, fit$obj$gr, lower=lower, upper=upper))
  res <- max(abs(sapply(fits,function(x)c((x$par-fit$opt$par)/fit$opt$par,nll=(x$objective-fit$opt$objective)/fit$opt$objective))))
  if(res>0.001){
    out<- t(sapply(fits,function(x)c(x$par,nll=x$objective, conv=x$conv)))
    options(width=200)
    cat(fit$label,"\n",file="jit.tab", append=TRUE)
    cat(sub("^[0-9]*","  ",capture.output(as.data.frame(out))), file = 'jit.tab', append=TRUE, sep = '\n')
  }
  res
}

dat <- read.table("Y.tab", head=FALSE)
mymap <- list(logSdObs=factor(rep(1,ncol(dat))))

resflag<-FALSE

pdf("res.pdf")
  mod <- list()
  mod[[length(mod)+1]] <- runit(mode=0, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod0-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1011, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi1-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1101, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi2-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1110, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi3-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1100, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi1-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1010, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi2-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=1001, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi3-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=2, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod2-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=3, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod3-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  mod[[length(mod)+1]] <- runit(mode=4, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod4-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  mod[[length(mod)+1]] <- runit(mode=4, trans=1/3, res=resflag,map=mymap, cut.data=10, label="Mod4-cubrt-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  mod[[length(mod)+1]] <- runit(mode=4, trans=1/2, res=resflag,map=mymap, cut.data=10, label="Mod4-sqrt-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  mod[[length(mod)+1]] <- runit(mode=5, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod5-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  mod[[length(mod)+1]] <- runit(mode=6, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod6-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
dev.off()



dat <- read.table("Y.tab", head=FALSE)
mymap <- list(logSdObs=factor(rep(1,ncol(dat))))


res <- as.data.frame(do.call(rbind, lapply(mod, function(m)c(m$label, round(m$logLik,2), round(m$AICc,2), round(m$AIC,2), m$conv))))

cv<-lapply(mod, function(m)cv.rmse(year=10, cv.scale=log, mode=m$call$mode, transCode=m$call$transCode, label=m$label, cut.data=m$call$cut.data, map=m$call$map,
                                   lowerLog=m$call$lowerLog, upperLog=m$call$upperLog, lowerLogit=m$call$lowerLogit, upperLogit=m$call$upperLogit))

jit <- sapply(mod, function(m)jitfun(m))

res<-cbind(res,do.call(rbind,cv),round(jit,3))

names(res)<-c("Label", "nlogLik", "AICc","AIC","Conv_all", "RMSE-CV", "Conv_rate_CV", "jit")

options(width=200)
cat(sub("^[0-9]*","  ",capture.output(res)), file = 'res.tab', sep = '\n')
options(width=80)
