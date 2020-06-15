library(stockassessment)

fit = fitfromweb("haddock_NS_2020_v002")

sw = fit$data$stockMeanWeight

write.table(sw,file="Y.tab",row.names=FALSE,col.names=FALSE)
                  
logN = t(fit$pl$logN)

write.table(exp(logN),file="N.tab",row.names=FALSE,col.names=FALSE)

