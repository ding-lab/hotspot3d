#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
y = read.table(args[1])
z = read.table(args[2])

RD<-y[[2]]
ID<-y[[1]]

x0<-z[[1]] 
y0<-z[[3]] 
x1<-z[[2]]+1 
y1<-z[[3]] 

pdf(args[3],width=23.6,height=13.3)
par(mar=c(8,5,5,1))
barplot(RD,names.arg=ID,ylab="Reachabilty Distance (A)",main=paste("Reachability Plot: Epsilon=",args[4],"MinPts=",args[5]),col="Red", border=NA, space=0, las=2, cex.names=0.4)
segments (x0,y0,x1,y1)
dev.off()

# args: 1-RD.out, 2-clusters.out, 3-pdf_file_name, 4-epsilon, 5-MinPts, 6-cutoff
