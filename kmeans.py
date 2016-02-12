__author__ = 'KATRINA'

#from example: http://mlpy.sourceforge.net/docs/3.3/cluster.html

import sys

import numpy as np
import matplotlib.pyplot as plt
import mlpy
#import sklearn
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D

#GET INPUT DATA
data_path =sys.argv[1]

d = np.loadtxt(data_path, delimiter="\t",skiprows=1)
dt = d.transpose()

estimator = KMeans(n_clusters=5, init='k-means++', n_init=10, max_iter=300, tol=0.0001, precompute_distances='auto', verbose=0, random_state=None, copy_x=True, n_jobs=1)

fignum = 1
fig = plt.figure(fignum, figsize=(4, 3))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
plt.cla()
estimator.fit(dt)
labels=estimator.labels_
ax.scatter(d[:, 3], d[:, 0], d[:, 2])#, c=labels.astype(np.float))

ax.w_xaxis.set_ticklabels([])
ax.w_yaxis.set_ticklabels([])
ax.w_zaxis.set_ticklabels([])
ax.set_xlabel('Petal width')
ax.set_ylabel('Sepal length')
ax.set_zlabel('Petal length')

plt.show()

'''
#COULD NOT FIND MLPY.KMEANS FUNCTION -- MAYBE REVISIT THIS IF SCIKIT METHOD DOESN"T WORK
cls, means, steps = mlpy.kmeans(d, k=7, plus=True)
#^^^requires that k is specified, which could be problematic for real data


fig = plt.figure(1)
plot1 = plt.scatter(d[:,0], x[:,1],c=cls, alpha=0.75)
plot2 = plt.scatter(means[:,0], means[:,1], c=np.unique(cls), s=128, marker='d')
fig.savefig(data_path+"outputfigure.png")
#plt.show()
'''