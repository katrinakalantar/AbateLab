__author__ = 'KATRINA'


import numpy as np
import sys
import fastcluster
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

if sys.argv[1] == "-h":
    print("usage: clusterData.py {input data matrix} {input refseq taxons} {output file name} ")

#GET INPUT DATA
data_path = sys.argv[1]
refseq_taxons_path = sys.argv[2]
output_file = sys.argv[3]

d = np.loadtxt(data_path, delimiter="\t",skiprows=1)
dt = d.transpose()

refseq_taxons = open(refseq_taxons_path,'r').readlines()

result3 = fastcluster.linkage_vector(dt, method='centroid', metric='euclidean')
np.savetxt(output_file+"_dendroMatrix.txt",result3,delimiter='\t')

fig1 = plt.figure(figsize=(10, 35))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(result3, color_threshold=1, show_leaf_counts=True, orientation='right')#labels=refseq_taxons,) #, truncate_mode='level',p=2 )
#plt.show()
fig1.savefig(output_file+"_dendro.png",dpi=900)