__author__ = 'KATRINA'
'''

viewAsGraph.py [mash distance file]
- use this to view the distance file as a graph - cannot run on abatelab computer bc of package dependencies
- only works for small sets of files when run locally

'''
import numpy as np
import pandas as pd
import sys
import networkx as nx


def calculateDistanceDF(mash_dist_file):
    print("inside calculateDistanceDF")
    pmF = open(mash_dist_file,'r').readlines()
    print(mash_dist_file)
    pairwise_dict = {}
    for line in pmF:
        l = line.split('\t')
        #distances are 1-(fraction of shared kmers in hash)
        #pairwise_dict[(l[0].split('/')[-1],l[1].split('/')[-1])] = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))
        pairwise_dict[(l[0].split('/')[-1],l[1].split('/')[-1])] = int((1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1])))<=.97)
    ser = pd.Series(list(pairwise_dict.values()),index=pd.MultiIndex.from_tuples(pairwise_dict.keys()))
    df = ser.unstack().fillna(0)
    return [df, pairwise_dict]


mash_distanceFile = sys.argv[1]
df, pairwise_dict = calculateDistanceDF(mash_distanceFile)

G = nx.from_numpy_matrix(np.matrix(df))
nx.draw(G,'graph.png')
