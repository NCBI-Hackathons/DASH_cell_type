"""Use parameters trained in grid search cross validation,
fit a model and compare prediction of test RNA-seq dataset."""

import sklearn
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import SGDClassifier
from collections import Counter


train = pd.read_table('final.txt')
del train['batch']
y = train['tissue']
del train['tissue']

mel = pd.read_table('GSE62526_Normalized_expression_values.txt.norm')
cll = pd.read_table('GSE66117_CLL_FPKM_values.txt.norm')
del mel['gene']
del cll['gene']

pca = PCA(40)
pca.fit(train)
X = pca.transform(train)

melpc = pca.transform(mel)
cllpc = pca.transform(cll)

#params from cross validation on training set
sgd = SGDClassifier('log', penalty='elasticnet', alpha=0.008, l1_ratio=0.4)
sgd.fit(X, y)

melpred = sgd.predict(melpc)
cllpred = sgd.predict(cllpc)

# all 29 are melanoma cancer
print melpred
#print Counter(melpred)
# Counter({'cancer_melanoma': 12, 'normal_ovary': 8, 'normal_lung': 7, 'normal_blood': 2})
# last 5 are controls, we got 1, acc 44.4%
print cllpred
#print Counter(cllpred)
# Counter({'cancer_leukemia': 48, 'normal_blood': 4}), acc 88.5%
