set.seed(123)
library(TMB)
# compile 
#compile("../../src/gmrf1.cpp")
dyn.load(dynlib("../../src/gmrf1"))

cat("###################\n",getwd(),"\n#####################")

runit <- function(mode=1, transCode=-1, res=FALSE, label=NULL, predict=0, cut=predict-1, silent=TRUE, cut.data=0, map=list(), lowerLog=-Inf, upperLog=Inf, lowerLogit=-Inf, upperLogit=Inf,  ...){
  cat("###################---\n",label,"\n#####################")
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
  Zorg <- as.matrix(read.table("Z.tab", head=FALSE))
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
  SSBorg<-rowSums(Yorg*Norg*Moorg)
  jac <- 0# -sum(log(abs(numDeriv:::grad(trans,Yorg[!is.na(Y)]))))
  
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

  addToMap<-function(x,map){
      for(i in 1:length(x)){
          map[[ names(x)[i] ]] <- x[[i]]
      }
      map
  }
    
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
    map<- addToMap(list(logPhi=factor(addmap)),map) 
    param$mu <- colMeans(data$Y, na.rm=TRUE)#numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    ran <- c("z")
  }
  if(mode==11){    
    param$logPhi <- c(0,0)
    param$mu <- colMeans(data$Y, na.rm=TRUE)#numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    ran <- c("z")
  }
  if(mode==13){    
    param$logPhi <- c(0,0)
    param$mu <- colMeans(data$Y, na.rm=TRUE)#numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$missing <- rep(0,sum(is.na(data$Y)))
    ran <- c("missing")
  }
  if(mode==15){    
    param$logPhi <- c(0,0)
    param$mu <- colMeans(data$Y, na.rm=TRUE)#numeric(ncol(data$Y))
    param$logSdProc <- rep(0,length(data$Y))
    map<- addToMap(list(logSdProc=factor(as.vector(col(data$Y)))),map) 
    param$missing <- rep(0,sum(is.na(data$Y)))
    ran <- c("missing")
  }
  #if(mode==16){
  #  data$logN<-log(Norg)
  #  data$Z<-Zorg  
  #  param$logPhi <- c(0,0)
  #  param$mu <- colMeans(data$Y, na.rm=TRUE)#numeric(ncol(data$Y))
  #  param$logSdProc <- 0
  #  param$logSdLogN <- c(0,0)
  #  param$missing <- rep(0,sum(is.na(data$Y)))
  #  ran <- c("missing")
  #}

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
    param$logPhi <- c(0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$logitRhoObs <- -1
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
  if(mode%in%c(4001,4010,4100,4110,4101,4011,4111)){
    data$mode <- 4
    idx<-as.numeric(strsplit(as.character(mode), "")[[1]])[-1]
    param$logitRho <- c(0,0,0)
    param$logitRho[idx==0]<- -10
    addmap<-1:3
    addmap[idx==0] <- NA
    map=addToMap(list(logitRho=factor(addmap)),map)

    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- c(0,0)
    
    param$logSdObs <- numeric(ncol(data$Y))
    param$omega <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$z <- rep(0,(nrow(Y)-1)+ncol(Y))
    if(idx[3]==0){ ## no cohort effect
        map <- addToMap(list(logSdProc=factor(c(1,NA)),z=factor(rep(NA,length(param$z)))),map )
        param$logSdProc[2] <- -10
    }
    param$logitRhoObs <- 0
    data$trans <- transCode
    stopifnot(transCode>=0)
    ran <- c("omega","z")
  }
  if(mode==10){    
    param$logPhi <- c(0,0,0)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$logLamBC <- 0
    ran <- c("z")
  }
  if(mode%in%c(10001,10010,10100,10110,10101,10011,10111)){
    data$mode <- 10
    idx<-as.numeric(strsplit(as.character(mode), "")[[1]])[-c(1:2)]
    param$logPhi <- c(0,0,0)
    param$logPhi[idx==0]<- -10
    addmap<-1:3
    addmap[idx==0] <- NA
    map<- addToMap(list(logPhi=factor(addmap)),map)
    param$mu <- numeric(ncol(data$Y))
    param$logSdProc <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$z <- matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))
    param$logLamBC <- 0
    ran <- c("z")
  }
  if(mode%in% c(120,121)){
    data$mode <- 12
    data$expIt <- mode %% 2  
    param$logitRho <- c(0,0)
    param$mu <- log(  pmax(diff(c(0,colMeans(data$Y,na.rm=TRUE))),1e-3) )  ##numeric(ncol(data$Y))
    ##if(data$expIt==0) param$mu <- exp(param$mu)
    
    param$logSdProc <- c(0)
    param$logSdObs <- numeric(ncol(data$Y))
    param$omega <- matrix(0,nrow=nrow(data$Y)+ncol(data$Y), ncol=ncol(data$Y))
    ran <- c("omega")
  }

  if(mode%in% c(140,141)){
    data$mode <- 14
    data$expIt <- mode %% 2
    data$N <- Norg
    param$logitRho <- c(0,0)
    param$mu <- log(  pmax(diff(c(0,colMeans(data$Y,na.rm=TRUE))),1e-3) )
    param$logSdProc <- c(0)
    param$logSdObs <- rep(-5,ncol(data$Y)); ##numeric(ncol(data$Y))
    param$omega <- matrix(0,nrow=nrow(data$Y)+ncol(data$Y), ncol=ncol(data$Y))
    param$lalpha <- rep(-5,ncol(data$Y)-1)
    ran <- c("omega")
  }
     
  if(mode==17){
    data$SSB <- SSBorg / mean(SSBorg)
    param$logWinf <- max(data$Y,na.rm=TRUE)
    param$logk <- log(0.1)
    param$logitRho0 <- 0
    param$logSdObs <- numeric(ncol(data$Y))
    param$logalpha <- -5
    param$logitRho1 <- c(10,10) ## Random walk
    param$logSdOmega <- c(-2,-2)
    param$omegak<- param$omegaw <- numeric(nrow(data$Y)+ncol(data$Y)) 
    ran <- c("omegaw","omegak")
    map<- addToMap(list(logitRho1=factor(c(NA,NA))),map)
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
    
    
  #opt <- nlminb(obj$par, obj$fn, obj$gr,lower=lower,upper=upper)


  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(eval.max=2000, iter.max=1000, rel.tol=1e-10),lower=lower,upper=upper)
  #for(i in seq_len(3)) { # Take a few extra newton steps 
  #  g <- as.numeric( obj$gr(opt$par) )
  #  h <- optimHess(opt$par, obj$fn, obj$gr)
  #  opt$par <- opt$par - solve(h, g)
  #  opt$objective <- obj$fn(opt$par)
  #}
  
  sdr <- sdreport(obj, opt$par)
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

cv2.rmse <- function(year=10, cv.scale=identity, ...){
    sq.error<-function(p){fit<-runit(...,predict=p, cut=p-2); c((tail(cv.scale(fit$ssbobs),1)-tail(cv.scale(fit$ssbpred),1))^2,as.integer(fit$conv))}
    ret<-rowMeans(Vectorize(sq.error)((1:year)+1))
    ret[1]<-sqrt(ret[1])
    ret
}

cv3.rmse <- function(year=10, cv.scale=identity, ...){
    sq.error<-function(p){fit<-runit(...,predict=p, cut=p-3); c((tail(cv.scale(fit$ssbobs),1)-tail(cv.scale(fit$ssbpred),1))^2,as.integer(fit$conv))}
    ret<-rowMeans(Vectorize(sq.error)((1:year)+2))
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
  mod[[length(mod)+1]] <- runit(mode=0, trans=0, res=resflag, map=mymap, cut.data=0, label="Mod0-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1011, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi1-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=11, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod11-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=11, trans=0, res=resflag, map=mymap, cut.data=0, label="Mod11-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=13, trans=0, res=resflag, cut.data=0, label="Mod13-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=15, trans=0, res=resflag, cut.data=0, label="Mod15-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1101, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi2-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1110, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1noPhi3-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1100, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi1-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1010, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi2-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=1001, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod1Phi3-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10, res=resflag, map=mymap, cut.data=10, label="Mod10-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10011, res=resflag, map=mymap, cut.data=10, label="Mod10noPhi1-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10101, res=resflag, map=mymap, cut.data=10, label="Mod10noPhi2-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10110, res=resflag, map=mymap, cut.data=10, label="Mod10noPhi3-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10100, res=resflag, map=mymap, cut.data=10, label="Mod10Phi1-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10010, res=resflag, map=mymap, cut.data=10, label="Mod10Phi2-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=10001, res=resflag, map=mymap, cut.data=10, label="Mod10Phi3-boxcox-constVar")
  #mod[[length(mod)+1]] <- runit(mode=2, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod2-log-constVar")
  #mod[[length(mod)+1]] <- runit(mode=3, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod3-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  #mod[[length(mod)+1]] <- runit(mode=4, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod4-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  #mod[[length(mod)+1]] <- runit(mode=4, trans=1/3, res=resflag,map=mymap, cut.data=10, label="Mod4-cubrt-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  #mod[[length(mod)+1]] <- runit(mode=4, trans=1/2, res=resflag,map=mymap, cut.data=10, label="Mod4-sqrt-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  #mod[[length(mod)+1]] <- runit(mode=6, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod6-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  #mod[[length(mod)+1]] <- runit(mode=4110, trans=0, res=resflag,map=mymap, cut.data=10, label="Mod4-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-4, upperLogit=4)
  ##mod[[length(mod)+1]] <- runit(mode=121, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod12.1-log-constVar")
  ##mod[[length(mod)+1]] <- runit(mode=120, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod12.0-log-constVar")
  mod[[length(mod)+1]] <- runit(mode=121, trans=0, res=resflag, map=mymap, cut.data=0, label="Mod12.1-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-5, upperLogit=5)
  mod[[length(mod)+1]] <- runit(mode=120, trans=0, res=resflag, map=mymap, cut.data=0, label="Mod12.0-log-constVar", lowerLog=-5, upperLog=5, lowerLogit=-5, upperLogit=5)
  ##mod[[length(mod)+1]] <- runit(mode=141, trans=0, res=resflag, map=list(logSdObs=factor(rep(NA,ncol(dat)))), cut.data=10, label="Mod14.1-log-constVar", lowerLog=-10, upperLog=5, lowerLogit=-5, upperLogit=5)
  ##mod[[length(mod)+1]] <- runit(mode=140, trans=0, res=resflag, map=list(logSdObs=factor(rep(NA,ncol(dat)))), cut.data=10, label="Mod14.0-log-constVar", lowerLog=-10, upperLog=5, lowerLogit=-5, upperLogit=5)
 #mod[[length(mod)+1]] <- runit(mode=140, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod14.0-log-constVar", lowerLog=-10, upperLog=5, lowerLogit=-5, upperLogit=5)
  ##mod[[length(mod)+1]] <- runit(mode=17, trans=0, res=resflag, map=mymap, cut.data=10, label="Mod17-log-constVar",lowerLog=-8)
dev.off()

res <- as.data.frame(do.call(rbind, lapply(mod, function(m)c(m$label, round(m$logLik,2), round(m$AICc,2), round(m$AIC,2), m$conv))))

cv<-lapply(mod, function(m)cv.rmse(year=10, cv.scale=log, mode=m$call$mode, transCode=m$call$transCode, label=m$label, cut.data=m$call$cut.data, map=m$call$map,
                                   lowerLog=m$call$lowerLog, upperLog=m$call$upperLog, lowerLogit=m$call$lowerLogit, upperLogit=m$call$upperLogit))

cv2<-lapply(mod, function(m)cv2.rmse(year=10, cv.scale=log, mode=m$call$mode, transCode=m$call$transCode, label=m$label, cut.data=m$call$cut.data, map=m$call$map,
                                   lowerLog=m$call$lowerLog, upperLog=m$call$upperLog, lowerLogit=m$call$lowerLogit, upperLogit=m$call$upperLogit))

cv3<-lapply(mod, function(m)cv3.rmse(year=10, cv.scale=log, mode=m$call$mode, transCode=m$call$transCode, label=m$label, cut.data=m$call$cut.data, map=m$call$map,
                                   lowerLog=m$call$lowerLog, upperLog=m$call$upperLog, lowerLogit=m$call$lowerLogit, upperLogit=m$call$upperLogit))

jit <- sapply(mod, function(m)jitfun(m))

res<-cbind(res,do.call(rbind,cv),do.call(rbind,cv2),do.call(rbind,cv3),round(jit,3))

names(res)<-c("Label", "nlogLik", "AICc","AIC","Conv_all", "RMSE-CV", "Conv_rate_CV", "RMSE-CV2", "Conv_rate_CV2", "RMSE-CV3", "Conv_rate_CV3", "jit")

options(width=200)
cat(sub("^[0-9]*","  ",capture.output(res)), file = 'res.tab', sep = '\n')
options(width=80)
