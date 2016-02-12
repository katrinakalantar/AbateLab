__author__ = 'KATRINA'
'''

pairwiseAnalysis.py [mash distance file] [output file name]
- create a cluster map (heatmap) of the data in the mash pairwise distance file
- used this to visualize the potential for a sample to cluster

'''

import sys
import pandas as pd
import seaborn as sns

pairwise_mashFile = sys.argv[1]
output = sys.argv[2]

#read mash distance file into a dataframe
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
#fig.savefig(output+"sns_heatmap_output.png")

#display the clustermap - containing heatmap and hierarchical clustering data in the margin
sns_plot2 = sns.clustermap(df, cbar=True, cmap="summer", method='ward', figsize=(25,25))
fig2 = sns_plot2.savefig(output+"sns_clustermap_output.png",dpi=1000)

