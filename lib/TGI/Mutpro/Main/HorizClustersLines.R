#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
y = read.table(args[1], sep = "\t")
z = read.table(args[2], sep = "\t")

RD<-y[[2]]
ID<-y[[1]]

z[z$V1==z$V2,"V3"] = 0.1 # show singletons at RD=0.1

y0<-z[[1]] 
x0<-z[[3]] 
y1<-z[[2]]+1 
x1<-z[[3]] 

Cluster<-z[[5]]

# adjust plot height according to the number of variants
newH = (24/300)*(length(ID)) # 300 is the number from MET where I used 23.6"
newH = round( newH, 1)
if (newH < 24){
  newH = 24
}

pdf(args[3],width=13.3,height=newH)
par(mar=c(8,5,5,1))
barplot(RD,names.arg=ID,main=paste("Reachability Plot: Epsilon=",args[4],"MinPts=",args[5]),col="Red", cex.names=0.4, horiz=TRUE,border=NA, space=0, las=2, xlab="Reachabilty Distance (A)") # ylab="Reachabilty Distance (A)"
segments (x0,y0,x1,y1)
text(x0,y1+0.5,Cluster, cex=0.4)
dev.off()

# args: 1-RD.out, 2-clusters.plot, 3-pdf_file_name, 4-epsilon, 5-MinPts