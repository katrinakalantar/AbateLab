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

#method to take cluster array and column headers and create dictionary of groups
def separateClusters(cluster_array, column_headers):
    max_cluster = max(cluster_array)
    print(max_cluster)
    #initialize a dictionary of groups
    groups = {}
    for i in range(1,max_cluster+1):
        groups[i]=[]
    #append each element to dictionary of groups
    for x in range(len(cluster_array)):
        groups[cluster_array[x]].append(column_headers[x].split('/')[-1])
    return groups

#input cluster dictionary, plot the composition of each group
def analyzeClusterComposition(cluster_dictionary):
    group_id = 0
    f = open(output+"_group_composition_dictionary",'w')


    #PLOT the distribution of species in each group

    fig = plt.figure()
    for c_list in cluster_dictionary:
        f2 = open(output +"_SAMPLECOMPOSITION_" + str(group_id),'w')
        print(str(c_list))
        group_id += 1
        x = cluster_dictionary[c_list]
        species_list = []
        #f2.write(str(group_id)+"-----\n")
        for i in x:
            species_list.append(i.split('-')[-2])
            f2.write(i+"\n")
        species_counts = Counter(species_list)
        f.write(str(species_counts))
        df=pd.DataFrame.from_dict(species_counts, orient='index')

        df.plot(kind='bar')
        plt.title('group ' + str(group_id)+' species distribution')
        locs,labels = plt.xticks()
        plt.setp(labels,rotation=0)
        plt.savefig(output+"graph"+str(group_id)+".png")



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

#CREATE the linkage file for clusters
#
#method = single, complete, centroid, average, weighted, median, ward
linkage_matrix = hierarchy.linkage(df, method='ward', metric='euclidean')
print('created linkage matrix')

dfm = df.as_matrix()
heatmapOrder=hierarchy.leaves_list(linkage_matrix)
orderedDataMatrix = dfm[heatmapOrder,:]
print(orderedDataMatrix)
orderedColHeaders = cols[heatmapOrder,:]

print(orderedColHeaders)
#t = 7, criterion='maxclust' -- setting the number of clusters to the known value
clusters=hierarchy.fcluster(linkage_matrix,19,'maxclust') ##THIS IS THE LINE I WANT TO USE WITH KNOWN CLUSTER NUMBERS
#clusters=hierarchy.fcluster(linkage_matrix,8,'inconsistent',depth=9)
print(clusters)

#print(clusters)
#print(cols)
cluster_dictionary = separateClusters(clusters, cols)
analyzeClusterComposition(cluster_dictionary)

with open(output+"_clusterDict_specificFilenames.txt",'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(list(cluster_dictionary.keys()))
    writer.writerows(zip(*cluster_dictionary.values()))

'''#commenting this out to hopefully save some time
# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    linkage_matrix,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
)
#plt.show()
plt.savefig(output + "_basic_dendrogram.png")
'''

#removed temporarily to speed up the distribution graphs
sns_plot3 = sns.clustermap(df, row_linkage = linkage_matrix, col_linkage = linkage_matrix, method='ward', cbar = True, cmap = "summer_r", figsize=(8,8))
#color options: ocean, summer,Spectral
plt.show()
fig3 = sns_plot3.savefig(output+"sns_clustermap_output.png")
