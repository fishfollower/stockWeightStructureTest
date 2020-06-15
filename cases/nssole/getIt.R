library(stockassessment)

fit = fitfromweb("Sole20_24_2020")

sw = fit$data$stockMeanWeight

write.table(sw[,-1],file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN[,-1]),file="N.tab",row.names=FALSE,col.names=FALSE)

##par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
