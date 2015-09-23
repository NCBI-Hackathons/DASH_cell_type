#PURPOSE: Combining Reference and NCI-60 data
#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library("GEOquery")
library(plyr)
library(ggplot2)

##REFERENCE##
REF="GDS596"
ref=getGEO(REF,destdir="~/Dash_Cell_Type")
deg=read.table("~/Dash_Cell_Type/deg.txt",header=FALSE)
colnames(deg)=c("IDENTIFIER")
data.matrix(Columns(ref)$tissue)

##subsetting based on nci-60 tissues
refeset <- GDS2eSet(ref,do.log2=F)
tissues=c("PB-CD4+T cells","PB-CD19+B cells","721-B-lymphoblasts","PB-CD56+NK cells",
          "PB-CD8+T cells","ovary","uterus","uterus corpus","skin",
          "hypothalamus","parietal lobe","prefrontal cortex","cerebellum",
          "kidney","lung","prostate")
tissues=tissues[-c(2,5,7,8)]
v=c()
for(x in tissues) {v=c(v,which(pData(refeset)[,2] %in% tissues))}
samples=pData(refeset)[unique(v),1]
b=c()
for(x in samples) {b=c(b,grep(x,names(Table(ref))))}
dim(Table(ref)[,c(1,2,b)])

## left-join the data ; keep only degs
ref_deg=join(deg,Table(ref)[,c(1,2,b)],by=c("IDENTIFIER"), type = "left", match = "first")
dim(ref_deg)

## romove nas
ref_deg_rmna = na.omit(ref_deg)
dim(ref_deg_rmna)


## add gene names as rownames
rownames(ref_deg_rmna)=ref_deg_rmna$IDENTIFIER

###TUMOR
gds=getGEO("GDS4296",destdir="~/Downloads")
deg=read.table("~/Downloads/deg.txt",header=FALSE)
colnames(deg)=c("IDENTIFIER")

## left-join the data ; keep only degs
exp_deg=join(deg,Table(gds),by=c("IDENTIFIER"), type = "left", match = "first")
dim(exp_deg)

## romove nas
exp_deg_rmna = na.omit(exp_deg)
dim(exp_deg_rmna)

## add gene names as rownames
rownames(exp_deg_rmna)=exp_deg_rmna$IDENTIFIER
###

#COMBING REFERENCE AND TUMOR
inter=intersect(rownames(exp_deg_rmna),rownames(ref_deg_rmna))
#scale NCI-60
exp_deg_rmna_mat=data.matrix(exp_deg_rmna)
exp_deg_rmna_mat_scale=scale(exp_deg_rmna_mat)
#scale reference
ref_deg_rmna_mat=data.matrix(ref_deg_rmna)
ref_deg_rmna_mat_scale=scale(ref_deg_rmna_mat)
df<-cbind(exp_deg_rmna_mat_scale[inter,],ref_deg_rmna_mat_scale[inter,-c(1,2)])
dim(df)
## transform and scale data
df_trans = t(data.matrix(df[,-(1:2)]))
dim(df_trans)
df_trans_mat = data.matrix(df_trans)
dim(df_trans_mat)
df_trans_mat_scale=df_trans_mat
        # scale(df_trans_mat)
dim(df_trans_mat_scale)

### batch adjust-still working on it
df_trans_mat_scale_trans=t(df_trans_mat_scale)
pheno=data.frame("sample"=rownames(df_trans_mat),"cancer"=groups,"batch"=c(rep(1,174),rep(2,24)),row.names=rownames(df_trans_mat))

mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)
mod
n.sv = num.sv(df_trans_mat_scale_trans,mod,method="leek")
n.sv
svobj = sva(df_trans_mat_scale_trans,mod,mod0,n.sv=2)
batch = pheno$batch
batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
###

svd = svd( df_trans_mat_scale )
U = svd$u
V = svd$v #PC loadings
D = svd$d
Z = df_trans_mat_scale%*%V #PCs
## get the top 4 pc
df_Top4_pc=U[,1:4]
colnames(df_Top4_pc)=c("PC1","PC2","PC3","PC4")
## plot
par(mfrow=c(1,1))
groups1<-paste0("cancer_",Columns(gds)$tissue)
groups2 <- paste0("normal_",pData(refeset)[unique(v),2])
t<-c("brain","brain","brain","brain","brain","brain","brain","brain","ovary","ovary","prostate","prostate","blood","blood","blood","blood","blood","blood","lung","lung","skin","skin","kidney","kidney")
groups2<-paste0("normal_",t)
groups<-c(as.character(groups1),groups2)
p <- ggplot(data.frame(df_Top4_pc[,1:3]), aes(PC1, PC3)) 
p + geom_point(size=5, aes(colour = groups)) +
        theme_bw() +theme(text = element_text(size=20)) #+
       #  scale_shape_manual(values=1:nlevels(groups))

#writing out files for machine learning
v=data.frame("sample"=rownames(df_trans_mat),
             "cancer"=groups,"batch"=c(rep(1,174),rep(2,24)))

write.table(v,row.names=rownames(df_trans_mat)),
            file="~/Dash_Cell_Type/temp.txt",col.names=F,
            row.names=F,quote=F,sep="\t")
final=cbind("tissue"=v[,2],"batch"=v[,3],df_Top4_pc[,1:3],
            data.frame(df_trans_mat_scale))
write.table(final,"~/Dash_Cell_Type/finalpc.txt",quote=F,row.names=F,col.names=T,sep="\t")

