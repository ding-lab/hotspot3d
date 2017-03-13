#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
d <- read.table(args[1], header=FALSE, sep = "\t") # data table with variant,RD,genomic_annotations,weight

segData = read.table(args[2], sep = "\t") # data table with start,stop,epsilon',clusterID

y0<-segData[[3]]
x0<-segData[[1]]+1
y1<-segData[[3]]
x1<-segData[[2]]+1
Cluster<-segData[[5]]

names(d)[1] <- "variant"
names(d)[2] <- "RD"
names(d)[10] <- "weight"

# adjust plot height according to the number of variants
newH = (24/300)*(length(d$variant)) # 300 is the number from MET where I used 23.6"
newH = round( newH, 1)
if (newH < 24){
  newH = 24
}
# replace RD=0 by 0.1 (so that we get a bar to color)
d$RD[d$RD==0] <- 0.1
d$weight <- log(d$weight)

library(ggplot2)

d$variant <- factor(d$variant, levels = d$variant) # avoid automatic sort

p <- ggplot(data=d, aes(x=variant, y=RD, fill=weight)) + geom_bar(stat="identity") + scale_fill_gradient(low="grey",high="red")

p <- p + coord_flip() + theme_bw() + theme(axis.text.y=element_text(size=6))  #+ theme_bw()

values <- c(1:nrow(segData))
for (i in values) {
  p <- p + geom_segment(x=x0[i], y=y0[i], xend=x1[i], yend=y1[i])
}

p <- p + annotate("text", x = x1+1, y = y0, label = Cluster, size = 2)

p <- p + ggtitle(paste("Reachability Plot with weights: Epsilon=",args[4],"MinPts=",args[5]))

ggsave(args[3], width = 13.3, height = newH, limitsize = FALSE)

# args: 1-RD.out, 2-clusters.plot, 3-pdf_file_name, 4-epsilon, 5-MinPts