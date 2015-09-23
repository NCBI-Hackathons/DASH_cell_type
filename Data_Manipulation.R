#http://genomicsclass.github.io/book/pages/pca_svd.html
setwd("~/Dash_Cell_Type")
install.packages("devtools")
library(devtools) 
install_github('dagdata','genomicsclass')
library(dagdata)
data(tissuesGeneExpression)
install.packages("rafalab")
data=df_trans_mat
pc<-prcomp(data)
names(pc)
group=final$tissue
plot(pc$x[, 1], pc$x[, 2], col = group, main = "PCA", xlab = "PC1", ylab = "PC2")
cx <- sweep(data, 2, colMeans(data), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = group, main = "SVD", xlab = "U1", ylab = "U2")
sv$v[1:5, 1:5]
pc$rotation[1:5, 1:5]
head(sv$d^2)
head(pc$sdev^2)
head(sv$d^2/(ncol(data) - 1))
plot(sv$d^2/sum(sv$d^2), xlim = c(0, 10), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
plot(sv$d^2/sum(sv$d^2), type = "b", pch = 16, xlab = "principal components", ylab = "variance explained")
svNoCenter <- svd(data)
plot(pc$x[, 1], pc$x[, 2], col = group, main = "PCA", xlab = "PC1", ylab = "PC2")
points(0, 0, pch = 3, cex = 4, lwd = 4)
plot(svNoCenter$u[, 1], svNoCenter$u[, 2], col = group, main = "SVD not centered", xlab = "U1", ylab = "U2")
#SVD on (genes vs samples) and (samples vs genes)
sv2 <- svd(t(data))
plot(sv2$u[, 1], sv2$u[, 2], col = group, main = "samples vs genes (typical PCA)", 
     xlab = "U1", ylab = "U2")
