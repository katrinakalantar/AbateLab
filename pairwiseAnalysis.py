__author__ = 'KATRINA'

import sys
#import numpy as np
#import matplotlib as plt
#import networkx
import pandas as pd
import seaborn as sns

pairwise_mashFile = sys.argv[1]
output = sys.argv[2]

pmF = open(pairwise_mashFile,'r').readlines()
pairwise_dict = {}
for line in pmF:
    l = line.split('\t')
    a = float(l[4].strip().split('/')[0])
    b = float(l[4].strip().split('/')[1])
    pairwise_dict[(l[0],l[1])] = (a/b)

ser = pd.Series(list(pairwise_dict.values()),
                  index=pd.MultiIndex.from_tuples(pairwise_dict.keys()))
df = ser.unstack().fillna(0)
print(df.shape)
#sns_plot = sns.heatmap(df, cbar=True, cmap="YlGnBu")
#fig = sns_plot.get_figure()
#fig.show()
#fig.savefig(output+"sns_heatmap_output.png")

sns_plot2 = sns.clustermap(df, cbar=True, cmap="summer", method='ward', figsize=(25,25))
fig2 = sns_plot2.savefig(output+"sns_clustermap_output.png",dpi=1000)
#fig2.savefig(output+"sns_clustermap_output.png")
