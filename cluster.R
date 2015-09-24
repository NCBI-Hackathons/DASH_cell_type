library(gplots)
install.packages("matrixStats")
library(matrixStats)
install.packages("RColorBrewer")
library(RColorBrewer)
table=read.table("~/Downloads/finalpc.txt",sep="\t",header=T)
cancer=table[table$batch==1,-c(1:5)]
ref=table[table$batch==2,-c(1:5)]
both=data.matrix(table[,-(1:5)])
cancer_trans=t(cancer)

topVarGenes <- head(order(-rowVars(cancer_trans)),100)
topVarGenes
mat <- as.matrix(cancer_trans[ topVarGenes, ])
mat <- mat - rowMeans(mat)
colors1 <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
colors2<-c(rep("black",3),"red","blue",rep("black",11))
sidecols <- colors2[ subset(table, batch==1)[,"tissue"] ]
heatmap.2(mat,trace="none",col=colors1,ColSideColors=sidecols)
# legend(x="bottomright", legend=levels(subset(table, batch==1)[,"tissue"])[1:9], col=colors2,pch=4)
