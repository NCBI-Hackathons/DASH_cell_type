"""Script to normalize test single cell RNA sequencing dataset
and output common subset genes as in the final training set.
"""

import pandas as pd
import sklearn
from sklearn.preprocessing import StandardScaler


info = {'GSE57982': {'filename': 'GSE57982_primaryFpkmMatrix.txt',
                     'idcol': 'geneSymbol', # column name of gene IDs
                     'rmcol': ['geneID']},  # columns to remove
        'GSE62526': {'filename': 'GSE62526_Normalized_expression_values.txt',
                     'idcol': 'Gene',
                     'rmcol': []},
        'GSE66117': {'filename': 'GSE66117_CLL_FPKM_values.txt',
                     'idcol': 'Gene',
                     'rmcol': ['Description']}
        }

with open('final.txt') as fp:
    # take the list of final deferentially expressed genes
    cols = fp.readline().split()
    cols.remove('tissue')
    cols.remove('batch')
    deg = pd.DataFrame(data=cols, columns=['gene'])


for gse, gseinfo in info.iteritems():
    pre_data = pd.read_table(gseinfo['filename'])

    # use left join to subset test data
    pre_data = pre_data.rename(columns={gseinfo['idcol']: 'gene'})
    pre_data = pre_data.drop_duplicates(subset='gene')
    pre_data = pd.merge(deg, pre_data, how="left", on='gene')
    pre_data.fillna(0, inplace=True)
    genes = pre_data['gene']
    del pre_data['gene']
    for col in gseinfo['rmcol']:
        del pre_data[col]

    scaler = StandardScaler()
    norm = scaler.fit_transform(pre_data.values)  # norm across genes for each sample
    tnorm = pd.DataFrame(data=norm.T, columns=genes)  # transpose to let sample on the row
    tnorm.insert(0, column='gene', value=pre_data.columns)

    output_fn = gseinfo['filename'] + '.norm'
    tnorm.to_csv(output_fn, sep='\t', index=False)
