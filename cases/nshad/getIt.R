library(stockassessment)

fit = fitfromweb("haddock_NS_2020_v002")

sw = fit$data$stockMeanWeight

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

matur = fit$data$propMat

write.table(matur,file="Mo.tab",row.names=FALSE,col.names=FALSE)
