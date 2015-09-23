##########
JHU Data Science Hackathon, Genomics Section, Cell Type Hacking Group
##########

##########
GOAL

>Build a model to classify tissue and tumor types from GEO data from NCBI 
##########

##########
IN PRACTICE

>Exploring how to download, handle and curate GEO data, specifically NCI-60 cancer cell lines and reference normal cell lines, and how to use machine learning to train a data set and test the model on single cell genomics data and see how tissue and tumor types classify.
##########

##########
BACKGROUND

PMID: 25434802

Any given tumor is a mixture of normal cells and tumor cells. Single cell genomic technology can capture molecular data from a single tumor cell, bypassing molecular characterization of heterogenous tumor mixtures. However, this technology is very new and not very high throughput, so it would be valuable to know if traditional cell line data can correctly classify data from single cells, which strengthens how accurate the cell lines characterize tissue and tumor samples.
 
##########

##########
DATA

NCI-60 Cancer cell line panel http://www.ncbi.nlm.nih.gov/gds/4296

Reference-Large scale analysis of the human transcriptome http://www.ncbi.nlm.nih.gov/gds/596

deg.txt - got all protein coding genes (424) that are differentially expressed on the protein level in NCI-60
- http://129.187.44.58:7070/NCI60/query/globalExpression

Single Cell Melanoma and Leukemia data
GSE66117_CLL_FPKM_values.txt RNA-seq data of chronic lymphocytic leukemia http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66117
GSE62526_Normalized_expression_values.txt RNA-seq data of 29 melanoma cell lines http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62526

##########

##########
PROGRAMMING

Downloaded and handled GEO data in R

Model Training and Testing in iPython using the scikit-learn package

##########

##########
FINAL THOUGHTS

Summary:

Gene expression data from cancer cell lines can differentiate samples from different cancer types. We investigated the classification accuracy of cancer and normal cell line data, using a testing set of homogenous tumor data from single-cell tumor gene expression data. In producing a model using the cancer and normal cell line data as a training set, we ran into challenges that are common with these types of analyses, such as batch effects integrating multiplatform data and overfitting of the model. 

https://github.com/NCBI-Hackathons/DASH_cell_type


##########
