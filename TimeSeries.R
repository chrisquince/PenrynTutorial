library(vegan)
library(ggplot2)
library(grid)
library(plyr)
rm(list=ls())
ClusterK <- read.csv("clustering_gt1000_covR.csv",header=TRUE,row.names=1)
ClusterK <- t(ClusterK)
ClusterP <- ClusterK/rowSums(ClusterK)
Meta <- read.table("sharon_mappingR.txt",sep='\t',header=TRUE)

rownames(Meta) <- Meta$Sample
colnames(ClusterP) <- gsub("^","Cluster",colnames(ClusterP))

scg <- read.table("clustering_gt1000_scg.tsv",header=TRUE,row.names=1)
rownames(scg) <- gsub("^","Cluster",rownames(scg))
ClusterP <- t(ClusterP)
ClusterP <- ClusterP[rownames(scg),]

ClusterP75 <- ClusterP[rowSums(scg == 1)/36. > 0.75,]

ClusterP75 <- t(ClusterP75)

ClusterMeta75 <- cbind.data.frame(ClusterP75,Meta)
ClusterMeta75$Sample <- NULL

ClusterMelt <- melt(ClusterMeta75)
colnames(ClusterMelt) <- c("Day","MAG","Freq")
p <- ggplot(data=ClusterMelt,aes(x=Day,y=Freq,colour=MAG,group=MAG)) + geom_point()

p <- p + geom_line() + theme_bw() + ylab("Relative frequency")

p <- p + theme(axis.title = element_text(size=12, face="bold")) 

pdf("TimeSeries.pdf")
plot(p)
dev.off()

