__author__ = 'KATRINA'

'''

tnf_barebones.py [input_directory]
- "input_directory" must contain all the .fastq files for which you intend to calculate TNF
- this program 1. calculates the TNF for all reads within a file and
               2. calculates the std deviation of all TNF vectors
               3. outputs the mean std deviation for each file
- used to determine whether std deviation in TNF can distinguish pure files from contaminated files

output:
tnf_barebones.log : general log file
tnf_barebonesOut.o : tab delimited file contining-
                        filename, and the average of all standard deviation values for the file

'''

import numpy as np
import sys
import timeit
import os
import glob
from scipy.spatial import distance
from matplotlib.mlab import PCA
import itertools

#return the vector of frequencies of each kmer in a single read
def getKmerFreqs(line):

    alphabet = 'ACGT'
    kmers = [''.join(i) for i in itertools.product(alphabet, repeat = 4)]
    kmerD = {} #initialize kmer dictionary
    for k in kmers:
        kmerD[k] = 0

    total_kmer_count = 0
    for i in range(len(line)-4):
         total_kmer_count += 1
         kmerD[line[i:i+4]] += 1

    tnf_dict = {}
    for k in kmerD:
        tnf_dict[k] = float(kmerD[k]/float(total_kmer_count))
    s = sorted(tnf_dict)
    tnf_result_vector = []
    for value in s:
        tnf_result_vector.append(tnf_dict[value])
    return tnf_result_vector

'''
def getTNFdist(line1_tnf,line2_tnf):
    dst = distance.euclidean(line1_tnf,line2_tnf)
    return dst


def runPCA(all_kmer_vectors_array):
    pca_dictionary = {}
    results = PCA(np.asarray(all_kmer_vectors_array))
    for i in range(len(all_kmer_vectors_array)):
        pca_dictionary[i] = results.Y[i][0:10] #arbitrarily setting 10 as the number of PCs to use
    return pca_dictionary
'''


input_directory = sys.argv[1]
os.chdir(input_directory)
array_of_filenames = glob.glob('./*.f*q')
out_file = open('tnf_barebones.log','w')
out_file2 = open('tnf_barebonesOut.o','w')

for input_file in array_of_filenames: #for each file in the input directory
    print(input_file)
    out_file.write(input_file + "\n")
    start = timeit.default_timer()
    f = open(input_file,'r').readlines()
    tnf_array = []
    for line_num in range(len(f)):
        tnf = []
        if line_num%4 == 1:
            tnf = getKmerFreqs(f[line_num]) #get the TNF for each read
            tnf_array.append(tnf)

    print(len(tnf_array))
    print(tnf_array)

    arr = np.asarray(tnf_array)
    std = np.std(arr, axis=0) #calculate the std dev for all TNF arrays in the file
    print(std)
    out_file.write(str(std))
    out_file.write('\n')
    out_file.write(str(np.mean(std))) #output the mean of this value
    out_file.write("\n\n")

    out_file2.write(input_file + '\t' + str(np.mean(std))+'\n')

    stop = timeit.default_timer()
    print("time to complete "+str(input_directory) +": " + str(stop - start)) #print time stats for each file



