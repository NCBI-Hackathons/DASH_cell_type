"""Script to do grid search for SVM and SGDClassifier."""

import sklearn
from sklearn import svm, grid_search
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest
from sklearn.linear_model import SGDClassifier
import pandas as pd 


train = pd.read_table('final.txt')
del train['batch']
y = train['tissue']
del train['tissue']
"""
y2 = []
for i in y:
    if i == 'cancer_melanoma':
        y2.append(i)
    elif i == 'normal_skin':
        y2.append(i)
    elif i == 'cancer_leukemia':
        y2.append(i)
    elif i == 'normal_blood':
        y2.append(i)
    elif 'normal' in i:
        y2.append('normal_other')
    else:
        y2.append('cancer_other')
#y = y2
"""
pca = PCA(40)
train = pca.fit_transform(train)

"""
# params grid for SVM
grid = [{'C': [1, 3, 5, 8, 10, 12, 15, 20, 50, 100, 500, 1000, 10000], 'kernel': ['linear', 'rbf'],
        'gamma': [0.0000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1]}
        ]
svr = svm.SVC()
"""

# params grid for SGD
grid = [{'loss' : ['log', 'hinge'], 'alpha': [0.005, 0.008, 0.01, 0.012, 0.05, 0.1, 0.15, 0.2, 0.5, 1, 10], 'penalty': ['l1', 'elasticnet'], 'l1_ratio': [0.2, 0.25, 0.3, 0.4, 0.45, 0.48, 0.5, 0.52,  0.55]}]
sgd = SGDClassifier()

#clf = grid_search.GridSearchCV(svr, grid, cv=3)
clf = grid_search.GridSearchCV(sgd, grid, cv=4)
clf.fit(train, y)
# Best Score:
# 0.848484848485
# Best Params:
# {'penalty': 'elasticnet', 'alpha': 0.008, 'loss': 'log', 'l1_ratio': 0.4}

print clf.grid_scores_
print 'Best Score:'
print clf.best_score_
print 'Best Params:'
print clf.best_params_
