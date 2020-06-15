library(stockassessment)

fit = fitfromweb("Herring-corrObs-update")

sw = fit$data$stockMeanWeight

bad = 1960:1982 - 1959

write.table(sw[-bad,],file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN[-bad,]),file="N.tab",row.names=FALSE,col.names=FALSE)

##par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
