library(stockassessment)

fit = fitfromweb("NEAcod-2020")

sw = fit$data$stockMeanWeight

good = !colnames(sw) %in% c(12:15) 

sw = sw[,good]

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)[,good]

matur = fit$data$propMat[,good]

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)

##par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
