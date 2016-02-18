__author__ = 'KATRINA'

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

'''

tnfvmash.py
- used to generate plots (scatter and density curve) comparing the mash dist to the tnf dist
- file names are hard coded within the script, this was for investigative purposes only
- the file names hard coded are from separated columns based on the file
    kmerTNF.o (output of early version of kmerTNFcombo.py)
    
* plots generated with this approach were less informative than expected,
  and therefore I stopped development in this direction, script may not be robust with other inputs

'''

#hard coded file names
kmer_file = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\TNF\\mashDist2"
tnf_file = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\TNF\\tnfDist2"
labels_File = "C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\TNF\\labels2"

kf = open(kmer_file,'r').readlines()
tf = open(tnf_file,'r').readlines()
lf = open(labels_File, 'r').readlines()

ka = [];ta = [];la = [];colors = []
ka_same = []; ka_dif = []; ta_same = []; ta_dif = []

for i in range(len(kf)):
    ka.append(kf[i])
    ta.append(float(tf[i]))
    la.append(lf[i])

    if lf[i].split(':')[0].split('-')[2] == lf[i].split(':')[1].split('-')[2]:
        colors.append('g')
        ka_same.append(kf[i])
        ta_same.append(float(tf[i]))
    else:
        colors.append('r')
        ka_dif.append(kf[i])
        ta_dif.append(float(tf[i]))

ta_max = max(ta)

ta_norms = []
ta_normd = []
for i in range(len(ta_same)):
    ta_norms.append(ta_same[i]/ta_max)
for i in range(len(ta_dif)):
    ta_normd.append(ta_dif[i]/ta_max)
'''
#CREATE SCATTER PLOTS
fig1 = plt.figure()
plt.scatter(ka_same, ta_norms,c='green')
plt.show()

fig2 = plt.figure()
plt.scatter(ka_dif, ta_normd, c= 'red')
plt.show()

plt.scatter(ka,ta_norm,c=colors,alpha=0.1)
#plt.scatter(ka,ta,c=colors,alpha=0.1)
'''


#CUMULATIVE DISTRIBUTION PLOT
dataA1 = ka_same; dataA2 = ta_norms; dataB1 = ka_dif
dataB2 = ta_normd
sorted_dataA1 = np.sort(dataA1)  # Or data.sort(), if data can be modified
sorted_dataA2 = np.sort(dataA2)
sorted_dataB1 = np.sort(dataB1)
sorted_dataB2 = np.sort(dataB2)
# Cumulative distributions:
fig3 = plt.figure()
plt.step(sorted_dataA1, np.arange(sorted_dataA1.size),c='o',label='kmer same sp')  # From 0 to the number of data points-1
plt.step(sorted_dataA2, np.arange(sorted_dataA2.size),c='b', label='tnf same sp')  # From the number of data points-1 to 0

plt.step(sorted_dataB1, np.arange(sorted_dataB1.size),c='g',label='kmer dif sp')  # From 0 to the number of data points-1
plt.step(sorted_dataB2, np.arange(sorted_dataB2.size),c='r', label='tnf dif sp')  # From the number of data points-1 to 0
plt.legend()

plt.title('density curve for same and different species')
ax3 = plt.axes()
ax3.set_yscale('log')
plt.show()
