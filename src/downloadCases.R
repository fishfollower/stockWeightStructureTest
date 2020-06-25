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

for( i in 1:length(stocks) ){

    cat(stocks[i],"...")
    
    thedir <- paste0("../cases/",names(stocks[i]))
    if(!dir.exists(thedir)){
        dir.create(thedir)
    }

    if(file.exists(paste0(thedir,"/Y.dat")) || overwrite ){
        fit <- fitfromweb(stocks[i],character.only=TRUE)
        sw <- fit$data$stockMeanWeight
        logN <- t(fit$pl$logN)
        matur <- fit$data$propMat
        logM <- fit$data$natMor
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

        sw <- sw[,goodAges]

        dsw <- apply(sw,2,diff)

        ## at most 1 repeated value pr year, AND always remove the last year (because it is often not observed, but set to avg. of last X years).
        
        goodYears <- c( apply(dsw,1,function(x) sum(x==0)<=1 ), FALSE )

        write.table(sw[goodYears,],file=paste0(thedir,"/Y.tab"),row.names=FALSE,col.names=FALSE)
        write.table(exp(logN)[goodYears,goodAges],file=paste0(thedir,"/N.tab"),row.names=FALSE,col.names=FALSE)
        write.table(matur[goodYears,goodAges],file=paste0(thedir,"/Mo.tab"),row.names=FALSE,col.names=FALSE)

        write.table(Z[goodYears,goodAges],file=paste0(thedir,"/Z.tab"),row.names=FALSE,col.names=FALSE)
        cat("done.\n")
        
    }   
}
