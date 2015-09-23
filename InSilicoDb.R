#Failed attempt to use this R package to download GEO data
#instead used GEOquery. Code for using GEOquery in GEO.R
InSilicoLogin("nick.giangreco@gmail.com","f381fec1724285b4c92cabf9c5c57575")
eset = getDatasetInfo("GSE32474", "GPL570")
print(eset)
scan = getDataset("GSE32474", "GPL570", norm = "SCAN", features = "gene")
print(exprs(scan)[1:10, 1:5])