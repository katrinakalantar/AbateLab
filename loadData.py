__author__ = 'KATRINA'

# January 8th 2015
# This script does the following:
# 1. Imports a file of mash-refseq hit vectors
# 2. Imports a file of refseq taxonomic classifications
# 3. Removes all rows containing zero values for all fastq files
# 4. Removes the corresponding refseq taxonomic classifications
# 5. Writes the new data matrix and refseq data to files

import numpy as np
import sys

if sys.argv[1] == "-h":
    print("usage: loadData.py {input data matrix} {input refseq taxons} {output file name} ")


#GET INPUT DATA
data_path =sys.argv[1]
refseq_taxons_path = sys.argv[2]

#1/2. read in data and headers
headers = open(data_path).readlines()[0].strip().split('\t')
#the way the .sh script generates the columns leaves and empty first column, so usecols will skip this
d = np.loadtxt(data_path, delimiter="\t",skiprows=1,usecols=range(1,int(sys.argv[4])))

#2. get an array of the taxons for reference
refseq_taxons = open(refseq_taxons_path,'r').readlines()

#print(d)
#print(refseq_taxons)

y = np.sum(d,axis=1) #sum of all rows in data table
rows_2be_removed = np.where(y==0)[0] #get index values of rows with now entries

#3/4. delete rows with only zeros (ie row sum == 0)
rows_already_deleted = 0 #keep track of rows already deleted to adjust indexing
for i in rows_2be_removed:
    d = np.delete(d, i-rows_already_deleted, axis=0)
    refseq_taxons = np.delete(refseq_taxons, i-rows_already_deleted, axis=0)
    rows_already_deleted += 1

#print(d)
#print(refseq_taxons)

#5. save the files
np.savetxt(sys.argv[3], d, delimiter='\t',header='\t'.join(headers))
open(sys.argv[3]+"_refseqTaxons_modified",'w').write(''.join(refseq_taxons))