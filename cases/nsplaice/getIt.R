library(stockassessment)

fit = fitfromweb("ple420_test")

sw = fit$data$stockMeanWeight[,1:10]

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

 par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
