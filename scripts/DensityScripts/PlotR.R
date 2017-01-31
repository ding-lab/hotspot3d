#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
y = read.table(args[1])

RD<-y[[2]]
ID<-y[[1]]

pdf(args[2],width=23.6,height=13.3)
par(mar=c(8,5,5,1))
barplot(RD,names.arg=ID,ylab="Reachabilty Distance (A)",main=paste("Reachability Plot: Epsilon=",args[3],"MinPts=",args[4]),col="Red", border=NA, space=0, las=2, cex.names=0.4)
dev.off()



