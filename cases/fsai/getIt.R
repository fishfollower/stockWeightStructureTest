library(stockassessment)

fit = fitfromweb("sam-tmb-fsaithe-2017-01")

sw = fit$data$stockMeanWeight

goodA = colnames(sw) %in% 3:14 

sw = sw[,goodA]

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)[,goodA]

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

matur = fit$data$propMat[,goodA]

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)

par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
