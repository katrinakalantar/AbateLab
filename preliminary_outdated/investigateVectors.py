__author__ = 'KATRINA'

'''
OUTDATED - this was used to understand the composition of vectors resulting from MASH distance comparison to RefSeq DB
this methods was abandoned in pursuit of the pairwise MASH distance approach
'''

import numpy as np
import sys
import matplotlib.pyplot as plt

# Create histograms for two different applications
# 1. plot the sum of hits for each vector loci in the refseq db -
#   this gives an idea of which loci are most commonly  hit - not sure how useful this is
# 2. For each genome, what was the maximum number of hits? plot the distribution of maximum hits

if sys.argv[1] == "-h":
    print("usage: investigateVectors.py {input data matrix} {input refseq taxons}")

#GET INPUT DATA
data_path =sys.argv[1]
refseq_taxons_path = sys.argv[2]

# read in data and headers
headers = open(data_path).readlines()[0].strip().split('\t')
d = np.loadtxt(data_path, delimiter="\t",skiprows=1)

# get an array of the taxons for reference
refseq_taxons = open(refseq_taxons_path,'r').readlines()

y = np.sum(d,axis=1) #sum of all rows in data table (ie hits at refseq loci across genomes)
x_axis_integers = []
for i in range(1,len(y)+1):
    x_axis_integers.append(i)
plot_dict = dict(zip(x_axis_integers,y))

#PLOT MAX values of hit for each genome
#PLOT Y TO MAKE THE SECOND HISTOGRAM (ACTUALLY A BAR CHART) THAT I WAS INTERESTED IN
fig, ax = plt.subplots()
ax.bar(plot_dict.keys(), plot_dict.values(),color='g')
fig.set_size_inches(18.5, 10.5)
ax.set_title('Sum of each row in vector')
fig.savefig(sys.argv[3]+"_img2.png")

#PLOT MAX values of hit for each genome
ymax = np.amax(d,axis=0)
print(ymax)
hist_dict = {}

for i in range(int(max(ymax))+1):  #create a dictionary of all possible values from 0 to max in a column (genome)
    hist_dict[i]=0
for j in np.array(ymax):  #count up each occurrence of a max value and add it to the dictionary
    hist_dict[j]=hist_dict[j]+1
print(hist_dict)

#PLOT hist_dict as a histogram to SATISFY THING #1 to do in my notebook
fig, ax = plt.subplots()
ax.bar(hist_dict.keys(), hist_dict.values(),color='g')
ax.set_title('Distribution of Max Hits in RefSeq DB')
fig.savefig(sys.argv[3]+"_img1.png")
