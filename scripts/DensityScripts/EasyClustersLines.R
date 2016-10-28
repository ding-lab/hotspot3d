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

Cluster<-z[[5]]
xpos<- c(0:(length(ID)-1))

pdf(paste("./Results/",args[3]),width=23.6,height=13.3)
par(mar=c(8,5,5,1))
barplot(RD,ylab="Reachabilty Distance (A)",main=paste("Reachability Plot: Epsilon=",args[4],"MinPts=",args[5]),col="Red", border=NA, xaxt="n", xlab="" , space=0)#names.arg=ID, space=0, las=2, cex.names=0.4
segments (x0,y0,x1,y1)
segments (xpos+0.5,0,xpos+0.5,-0.05)
segments (xpos+0.5,-0.05,xpos+0.2,-0.08)
text(x1+1,y0,Cluster, cex=0.4)

#op<- par( mar = c(10,1,0,1))
text( x=xpos-1.35, y=-0.22 ,labels=ID, srt=45, cex=0.3, xpd=TRUE) # par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
#par(op)
dev.off()

# args: 1-RD.out, 2-clusters.out, 3-pdf_file_name, 4-epsilon, 5-MinPts
