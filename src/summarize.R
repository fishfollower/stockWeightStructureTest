cases = list.dirs("../cases")
cases = setdiff(cases,c("../cases","../cases/sim1"))

restabs = lapply(cases,function(x) read.table(paste0(x,"/res.tab"),header=TRUE) )

cases.short = gsub("../cases/","",cases)
names(restabs) <- cases.short

numcols = rep(NA,ncol(restabs[[1]]))
for(i in 1:length(numcols)) numcols[i] = is.numeric(restabs[[1]][,i])

ares <- array(NA,dim= c(length(restabs),nrow(restabs[[1]]),ncol(restabs[[1]])-sum(!numcols)))
for(i in 1:length(restabs)) ares[i,,] <- as.matrix(restabs[[i]][,numcols])

sumtab = apply(ares,c(2,3),FUN=sum)
meantab = apply(ares,c(2,3),FUN=mean)

colnames(sumtab) <- colnames(meantab) <- colnames(restabs[[1]])[numcols]
rownames(sumtab) <- rownames(meantab) <- restabs[[1]]$Label

cat("meantab:\n")
print(meantab)

write.table(sumtab,"../sumtab.txt")
write.table(meantab,"../meantab.txt")
