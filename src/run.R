library(TMB)

runit <- function(mode=1, transCode=-1, res=FALSE, label=paste0(mode,",",deparse(substitute(trans))), ...){

  ## trans codes: -1 = identity, 0 = log, >0 = power [ x^(transCode) ]
  trans <- identity
  invtrans <- identity  
  if( transCode == 0 ){
      trans <- log
      invtrans <- exp
  } else if( transCode > 0){
      trans <- function(x) x^(transCode)
      invtrans <- function(x) x^(1/transCode)
  }
    
  # setup data 
  data <- list()
  cat("###################\n",getwd(),"\n#####################")
  Y <- as.matrix(read.table("Y.tab", head=FALSE))
  Y[abs(Y)<1.0e-12] <- NA
  jac <- -sum(log(abs(numDeriv:::grad(trans,Y[!is.na(Y)]))))
  Y <- trans(Y)  
    
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
  
  # compile 
  compile("../../src/gmrf1.cpp")
  dyn.load(dynlib("../../src/gmrf1"))

  # run model 
  obj <- MakeADFun(data,param,random=ran, DLL="gmrf1",...)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr <- sdreport(obj)
  pred <- as.list(sdr, report=TRUE, what="Est")$pred
  predSd <- as.list(sdr, report=TRUE, what="Std")$pred


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
    
  return(list(logLik=opt$objective+jac, obj=obj, residual=residual, opt=opt, sdr=sdr))
}

dat <- read.table("Y.tab", head=FALSE)
mymap <- list(logSdObs=factor(rep(1,ncol(dat))))

pdf("res.pdf")
  runit(mode=1, res=TRUE, map=mymap)
  runit(mode=1, trans=0, res=TRUE, map=mymap)  
  runit(mode=2, res=TRUE, map=mymap)
  runit(mode=2, trans=0, res=TRUE, map=mymap)    
  runit(mode=3, res=TRUE, map=mymap)
  runit(mode=3, trans=0, res=TRUE, map=mymap)
  runit(mode=4, trans=0, res=TRUE,map=mymap)
  runit(mode=4, trans=1/3, res=TRUE,map=mymap)
  runit(mode=4, trans=1/2, res=TRUE,map=mymap)
dev.off()

