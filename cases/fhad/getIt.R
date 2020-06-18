library(stockassessment)

fit = fitfromweb("NWWG2018_Had_addedsurveypoint")

sw = fit$data$stockMeanWeight

good = !rownames(sw) %in% 1957:1976
goodA = -1

sw = sw[good,goodA]

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)[good,goodA]

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

matur = fit$data$propMat[good,goodA]

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)

par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
