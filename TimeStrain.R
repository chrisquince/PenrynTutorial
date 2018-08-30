library(vegan)
library(ggplot2)
library(grid)
library(plyr)
library(reshape2)
library(getopt)
rm(list=ls())

spec = matrix(c('gammafile','g',1,"character",'metafile','m',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)

Gamma <- read.csv(opt$gammafile,header=TRUE,row.names=1)
GammaK <- Gamma
GammaP <- GammaK/rowSums(GammaK)
Meta <- read.table(opt$metafile,sep='\t',header=TRUE)

rownames(Meta) <- Meta$Sample
#colnames(GammaP) <- gsub("^","H",colnames(GammaP))
Meta <- Meta[rownames(GammaP),]

GammaMeta <- cbind.data.frame(GammaP,Meta)

GammaMelt <- melt(GammaMeta)
#[1] "Day"      "Sample"   "variable" "value"   
colnames(GammaMelt) <- c("Day","Sample","Strain","Freq")
print(colnames(GammaMelt))
p <- ggplot(data=GammaMelt,aes(x=Day,y=Freq,colour=Strain,group=Strain)) + geom_point()

p <- p + geom_line() + theme_bw() + ylab("Relative frequency")

p <- p + theme(axis.title = element_text(size=12, face="bold")) 

pdf("TimeSeries.pdf")
plot(p)
dev.off()

