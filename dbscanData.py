__author__ = 'KATRINA'


import sys
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import Counter
import matplotlib.pyplot as plt
import csv
from sklearn.cluster import DBSCAN


#input from command line
pairwise_mashFile = sys.argv[1]
output = sys.argv[2]

print(pairwise_mashFile)
print(output)

#CONVERT pairwise MASH distance .txt file into a dataframe representing the distance matrix
pmF = open(pairwise_mashFile,'r').readlines()
pairwise_dict = {}
for line in pmF:
    l = line.split('\t')
    a = float(l[4].strip().split('/')[0])
    b = float(l[4].strip().split('/')[1])
    pairwise_dict[(l[0],l[1])] = 1-(a/b)   #UPDATED THIS TO BE 1-(shared values) = dist - so the diagonals are 0
ser = pd.Series(list(pairwise_dict.values()),
                  index=pd.MultiIndex.from_tuples(pairwise_dict.keys()))
df = ser.unstack().fillna(0)

#SAVE DATAFRAME and columns, in order, for reference
cols = list(df.columns.values)
print(cols)
np.savetxt(output+'_dataframe.txt', df, delimiter="\t")
f = open(output+'_columnheaders.txt','w')
for c in cols:
    f.write(c)
df.to_csv(output+'_fullDF.txt',sep='\t',header=True)

db = DBSCAN(eps=0.5, min_samples=10,metric='precomputed').fit_predict(df)
#print(db.labels)
print(db)
for x in db:
    print(x)