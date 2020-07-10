library(stockassessment)

overwrite <- TRUE

stocks <- c("bw"="BW_2018",
            "fsai"="sam-tmb-fsaithe-2017-01",
            "neacod"="NEAcod-2020",
            "neasai"="NEA_sei_20v1",
            "nshad"="haddock_NS_2020_v002",
            "nsplaice"="nsplaice_test",
            "nssole"="Sole20_24_2020",
            "fhad"="NWWG2018_Had_addedsurveypoint",
            "mack"="MackWGWIDE2019v02",
            "neahad"="Nea_haddock_2019",
            "nscod"="nscod19_ass01_October",
            "nsher"="Herring-corrObs-update",
            "nssai"="NS_saithe_2018_rerun",
            "nswhit"="NSwhiting_2020_new_method_new1")

par(mfrow=n2mfrow(length(stocks)),mar=c(1,1,1,1))

for( i in 1:length(stocks) ){

    cat(stocks[i],"...")
    
    thedir <- paste0("../cases/",names(stocks[i]))
    if(!dir.exists(thedir)){
        dir.create(thedir)
    }

    if(file.exists(paste0(thedir,"/Y.dat")) || overwrite ){
        fit <- fitfromweb(stocks[i],character.only=TRUE)
        sw <- fit$data$stockMeanWeight
        orig.ages <- colnames(sw)
        logN <- t(fit$pl$logN)
        matur <- fit$data$propMat
        logM <- log(fit$data$natMor)
        if(all(fit$conf$keyLogFsta>=0)){
            logF <- t(fit$pl$logF)[,fit$conf$keyLogFsta+1]
        } else { ## not all F's are used => set to zero
            logF <- logM
            logF[] <- 0
            posF <- fit$conf$keyLogFsta > 0
            logF[,fit$conf$keyLogFsta[posF]+1] <- t(fit$pl$logF)[,fit$conf$keyLogFsta[posF]+1]
        }
        
        Z <- exp(logF) + exp(logM)

        ## remove ages with zero/negative values
        goodAges <- apply(sw,2,function(x) !any(x<=0))

        sw <- sw[,colnames(sw) %in% names(goodAges[goodAges])]
        
        dsw <- apply(sw,2,diff)
        
        ## at most 1 repeated value pr year, AND always remove the last year (because it is often not observed, but set to avg. of last X years).
        
        goodYears <- c( apply(dsw,1,function(x) sum(x==0)<=1 ), FALSE )

        sw <- sw[goodYears,]
        
        ## remove ages with repeated (imputed) values over time
        dsw <- apply(sw,2,diff)
        goodAges2 <- apply(dsw,2,function(x) sum(x==0)<=1 ) 
        sw <- sw[,colnames(sw) %in% names(goodAges2[goodAges2])]

        ## remove repeated age groups
        dsw2 <- apply(sw,1,diff)
        badAges <- apply(dsw2,1,function(x) sum(x==0)>1 )
        sw <- sw[,!colnames(sw) %in% names(badAges[badAges])]

        goodAges <- which( colnames(sw) %in% orig.ages )
        
        includesPlusGroup <- fit$conf$maxAge == as.numeric(max(colnames(sw)))

        info <- data.frame(includesPlusGroup = includesPlusGroup, minAge = min(as.numeric(goodAges)),maxAge = max(as.numeric(goodAges)))
        write.table(info,file=paste0(thedir,"/info.tab"))                                                                       
        
        write.table(sw,file=paste0(thedir,"/Y.tab"),row.names=FALSE,col.names=FALSE)
        write.table(exp(logN)[goodYears,goodAges],file=paste0(thedir,"/N.tab"),row.names=FALSE,col.names=FALSE)
        write.table(matur[goodYears,goodAges],file=paste0(thedir,"/Mo.tab"),row.names=FALSE,col.names=FALSE)

        write.table(Z[goodYears,goodAges],file=paste0(thedir,"/Z.tab"),row.names=FALSE,col.names=FALSE)

        onemat = matrix(1,nrow(Z),ncol(Z))
        Cpred <- exp(logF)/Z*(onemat-exp(-Z))*exp(logN)
        write.table(Cpred[goodYears,goodAges],file=paste0(thedir,"/C.tab"),row.names=FALSE,col.names=FALSE)
        
        cat("done.\n")
        matplot(sw,main=stocks[i],axes=FALSE)

        
    }   
}
