#PURPOSE: Playing with GEOquery
#Good tutorial from Sean Davis http://watson.nci.nih.gov/~sdavis/tutorials/publicdatatutorial/
#install.packages("GEOquery")
library("GEOquery")
#############
GDS="GDS4296" #NCI-60 Cancer cell line panel
REF="GDS596" #Large scale analysis of the human transcriptome
REFGSE="GSE1133"
GSE="GSE32474"
path="~/Dash_Cell_Type"
gds=getGEO(GDS,destdir=path)
gdseset=GDS2eSet(gds,do.log2=T)
gse=getGEO(GSE,destdir=path,GSEMatrix=F)
gsmlist=GSMList(gse)
ref=getGEO(REF,destdir=path) 
refgse=getGEO(REFGSE,destdir=path,GSEMatrix=F)
refgsmlist=GSMList(refgse)
REFGPL=Meta(ref)$platform
refeset=GDS2eSet(ref,do.log2=T)
dim(refeset)
names(Meta(gds)) #it's a list of metadata
dim(Table(gds)) #Table of values for samples at genes/ids
dim(Table(gds)[,grep("GSM",colnames(Table(gds)))]) #only columns with samples
names(Columns(gds)) #names of columns of GDS
num.samples=length(unique(Columns(gds)[,1]))
num.cell.line=length(unique(Columns(gds)[,2]))
num.tissue=length(unique(Columns(gds)[,3]))
num.disease.state=length(unique(Columns(gds)[,4]))
GPL=Meta(gds)$platform 
ma=GDS2MA(gds,do.log2=T)
save.image("~/Dash_Cell_Type/GEO.RData")
