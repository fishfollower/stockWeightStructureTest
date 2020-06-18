library(stockassessment)

fit = fitfromweb("nscod19_ass01_October")

sw = fit$data$stockMeanWeight

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

matur = fit$data$propMat

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)

par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
