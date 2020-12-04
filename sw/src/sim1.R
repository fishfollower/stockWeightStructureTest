set.seed(123)
nY<-60
nA<-10
n<-nA*nY

Y<-matrix(NA,nrow=nY, ncol=nA)
r<-as.vector(row(Y))
c<-as.vector(col(Y))

W.r<-matrix(0,nrow=n, ncol=n)
W.c<-matrix(0,nrow=n, ncol=n)
W.d<-matrix(0,nrow=n, ncol=n)

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

phi.r<-1
phi.c<-10
phi.d<-100
Q<-(diag(n)-phi.r*W.r-phi.c*W.c-phi.d*W.d)
S<-solve(Q)
x<-MASS:::mvrnorm(1,rep(0,n),Sigma=S)

Y<-exp(matrix(x+rnorm(length(x),0,sd=.1), ncol=nA))

onemat<-matrix(1,nrow=nrow(Y),ncol=ncol(Y))

write.table(Y, row.names=FALSE, col.names=FALSE, sep="\t", file="cases/sim1/Y.tab")
write.table(onemat, row.names=FALSE, col.names=FALSE, sep="\t", file="cases/sim1/N.tab")
write.table(onemat, row.names=FALSE, col.names=FALSE, sep="\t", file="cases/sim1/Mo.tab")

