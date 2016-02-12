__author__ = 'KATRINA'

import numpy as np
import pandas as pd
import sys
from collections import Counter
import timeit
import os
import glob
import subprocess
from scipy.spatial import distance
from matplotlib.mlab import PCA
import itertools

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


def getTNFdist(line1_tnf,line2_tnf):
    dst = distance.euclidean(line1_tnf,line2_tnf)
    return dst


def runPCA(all_kmer_vectors_array):
    pca_dictionary = {}
    results = PCA(np.asarray(all_kmer_vectors_array))

    for i in range(len(all_kmer_vectors_array)):
        pca_dictionary[i] = results.Y[i][0:10] #arbitrarily setting 10 as the number of PCs to use

    return pca_dictionary



input_directory = sys.argv[1]
os.chdir(input_directory)
array_of_filenames = glob.glob('./*.f*')
out_file = open('tnf_barebons.o','w')

for input_file in array_of_filenames:
    print('!!!')
    start = timeit.default_timer()
    #input_file = sys.argv[1]
    f = open(input_file,'r').readlines()
    tnf_array = []
    for line_num in range(len(f)):
        tnf = []
        if line_num%4 == 1:
            tnf = getKmerFreqs(f[line_num])
            #print(tnf)
            #print('---')
            tnf_array.append(tnf)
        #print(tnf_array)
    #print(tnf_array)
    print(len(tnf_array))



    #pca_result = runPCA(tnf_array)


    distances = np.zeros((len(tnf_array),len(tnf_array)))
    for i in range(len(tnf_array)):
        for j in range(len(tnf_array)):
            #print(tnf_array[i])
            #print(tnf_array[j])
            distances[i,j] = getTNFdist(tnf_array[i],tnf_array[j])
            #distances[i,j] = getTNFdist(pca_result[i],pca_result[j])
    print(distances)
    for i in distances:
        out_file.write(str(i))
        out_file.write('\n')
    out_file.write(str(distances))
    out_file.write('\n')

    stop = timeit.default_timer()
    print("time to complete "+str(input_directory) +": " + str(stop - start))



