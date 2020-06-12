library(TMB)

# setup data 
data <- list()

Y <- as.matrix(read.table("Y.tab", head=FALSE))
r <- as.vector(row(Y))
c <- as.vector(col(Y))
n <- length(r)
W.r <- W.c <- W.d <- matrix(0,nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    W.r[i,j]<- (r[i]==r[j])&(abs(c[i]-c[j])==1)
    W.c[i,j]<- (c[i]==c[j])&(abs(r[i]-r[j])==1)
    W.d[i,j]<- (((r[i]-r[j])==1)&((c[i]-c[j])==1))|(((r[i]-r[j])==(-1))&((c[i]-c[j]) ==(-1)))
  }      
}
diag(W.r)<- -rowSums(W.r)
diag(W.c)<- -rowSums(W.c)
diag(W.d)<- -rowSums(W.d)

data$Wr = W.r
data$Wc = W.c
data$Wd = W.d
data$Y = Y

# setup parameters 
param<-list()
param$logPhi<-c(0,0,0)
param$mu<-numeric(ncol(data$Y))
param$logSdProc<-0
param$logSdObs<-0
param$z<-matrix(0,nrow=nrow(data$Y), ncol=ncol(data$Y))

# compile 
compile("../../src/gmrf1.cpp")
dyn.load(dynlib("../../src/gmrf1"))

# run model 
obj<-MakeADFun(data,param,random="z", DLL="gmrf1")
opt<- nlminb(obj$par,obj$fn, obj$gr)
sdr<-sdreport(obj)
pred<-as.list(sdr, report=TRUE, what="Est")$pred
predSd<-as.list(sdr, report=TRUE, what="Std")$pred

pdf("res.pdf")
matplot(pred, type="l", ylim=range(data$Y))
#matplot(pred-2*predSd, , type="l", add=TRUE, lty="dotted")
#matplot(pred+2*predSd, type="l", add=TRUE, lty="dotted")
matplot(data$Y, add=TRUE)
dev.off()
