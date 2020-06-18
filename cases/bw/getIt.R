library(stockassessment)

fit = fitfromweb("BW_2018")

sw = fit$data$stockMeanWeight

good = 1:(nrow(sw)-1)

sw = sw[good,]

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)[good,]

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

matur = fit$data$propMat[good,]

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)

par(mfrow=n2mfrow(ncol(sw))); for(i in 1:ncol(sw)) plot(logN[,i],log(sw[,i]),main=i)
