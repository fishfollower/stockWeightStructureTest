library(stockassessment)

fit = fitfromweb("NSwhiting_2020_new_method_new1")

sw = fit$data$stockMeanWeight

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

## par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
