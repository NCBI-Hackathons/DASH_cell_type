#PURPOSE: Load and plot NCI-60 data
#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library("GEOquery")
path="~/Dash_Cell_Type"
install.packages(paste0(path,"/plyr"), type="source", repos=NULL)
library(ggplot2)


gds=getGEO("GDS4296",destdir=path)
unique(Columns(gds)$tissue)
unique(Columns(gds)$disease.state)
deg=read.table(paste0(path,"/deg.txt"),header=FALSE)
colnames(deg)=c("IDENTIFIER")
dim(deg)

## left-join the data ; keep only degs
#subsetting-going from 54675 features to 424 features (genes)
exp_deg=join(deg,Table(gds),by=c("IDENTIFIER"), type = "left", match = "first")
dim(exp_deg)

## romove nas
exp_deg_rmna = na.omit(exp_deg)
dim(exp_deg_rmna)

## add gene names as rownames
rownames(exp_deg_rmna)=exp_deg_rmna$IDENTIFIER

## transform and scale data
exp_deg_rmna_trans = t(data.matrix(exp_deg_rmna[,-(1:2)]))
exp_deg_rmna_trans_mat_scale=scale(exp_deg_rmna_trans)

## pca
svd = svd( exp_deg_rmna_trans_mat_scale )
U = svd$u
V = svd$v #PC loadings
D = svd$d
Z = exp_deg_rmna_trans_mat_scale%*%V #PCs

## get the top 4 pc
exp_Top4_pc=U[,1:4]
colnames(exp_Top4_pc)=c("PC1","PC2","PC3","PC4")

## plot
par(mfrow=c(1,1))
groups <- factor(Columns(gds)$tissue)
p <- ggplot(data.frame(exp_Top4_pc[,1:2]), aes(PC1, PC2)) 
p + geom_point(size=5, aes(colour = groups)) + theme_bw() +theme(text = element_text(size=20))

groups_disease <- factor(Columns(gds)$disease.state)
q <- ggplot(data.frame(exp_Top4_pc[,1:2]), aes(PC1, PC2)) 
q + geom_point(size=5, aes(colour = groups_disease)) + theme_bw() +theme(text = element_text(size=20))
