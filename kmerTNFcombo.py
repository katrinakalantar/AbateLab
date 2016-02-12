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


##NOTE: THIS IS STILL UNDER DEVELOPMENT - HAVE ONLY TESTED THE METRICS OBTAINED THROUGH ~/sandbox

def runMashSketch(file_string,k,s):
    print("inside runMashSketch")
    #print(file_string)
    cmd = mash+' sketch -o reference -u -k '+str(k)+' -s '+str(s)+' '+file_string
    print(cmd)
    os.system(cmd)
    #process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    #process = subprocess.Popen([mash,"sketch" "-o reference","-u","-k "+str(k), "-s "+str(s),file_string], shell=True, stdout=subprocess.PIPE)
    #process.wait()
    #print(process.returncode)
    print('finished runMashSketch')
    return

def runMashDist(file_string,k,s):
    print("inside runMashDist")
    cmd = mash+' dist reference.msh -s ' +str(s)+' '+file_string+'> '+'mash_pdist_'+str(k)+'_'+str(s)
    print(cmd)
    #os.system(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)
    print('finished runMashDist')
    return

def getKmerFreqs(input_file):
    f = open(input_file,'r').readlines()
    alphabet = 'ACGT'
    kmers = [''.join(i) for i in itertools.product(alphabet, repeat = 4)]
    kmerD = {} #initialize kmer dictionary
    for k in kmers:
        kmerD[k] = 0

    total_kmer_count = 0
    for line_num in range(len(f)):
        if line_num%4 == 1:
            line = f[line_num]
            #print(line)
            for i in range(len(line)-4):
                total_kmer_count += 1
                kmerD[line[i:i+4]] += 1

    #print(str(kmerD))
    new_kmerD = {}
    for k in kmerD:
        new_kmerD[k] = float(kmerD[k]/float(total_kmer_count))

    s = sorted(new_kmerD)
    kmer_result_vector = []
    for value in s:
        kmer_result_vector.append(new_kmerD[value])
    return kmer_result_vector

def calculateTNF(array_of_filenames):
    all_kmer_vectors = {}; all_kmer_vectors_array = []; labels = []
    for input_file in array_of_filenames:
        labels.append(str(input_file))#.split('-')[2])
        v = getKmerFreqs(input_file)
        all_kmer_vectors_array.append(np.asarray(v))
        all_kmer_vectors[input_file] = v
    return all_kmer_vectors_array, labels#all_kmer_vectors

def runPCA(all_kmer_vectors_array, labels):
    pca_dictionary = {}
    results = PCA(np.asarray(all_kmer_vectors_array))

    for i in range(len(labels)):
        pca_dictionary[labels[i]] = results.Y[i][0:10] #arbitrarily setting 10 as the number of PCs to use

    return pca_dictionary


def getTNFdist(pca_dictionary, file_name1, file_name2):
    #print(pca_dictionary.keys()[0:10])
    #print(file_name1)
    #print(file_name2)
    v1 = pca_dictionary[file_name1]
    #print(v1)
    v2 = pca_dictionary[file_name2]
    #print(v2)
    dst = distance.euclidean(v1,v2)
    return dst

def getClusters(mash_distanceFile, threshold, pca_dictionary):
    output_metrics = open("metrics.o",'w')
    unclustered = open("unclustered.o",'w')
    distFile = open(mash_distanceFile,'r')


    result_dictionary = {}
    temp_dictionary = {}
    all_file_names = {}

    line_count = 0
    number_of_groups = 0

    for line in distFile:
        line_count += 1

        l = line.split('\t')

        c1 = l[0].split('/')[-1]
        c2 = l[1].split('/')[-1]
        dist = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))

        #tnf_dist = getTNFdist(pca_tnf_vectors,l[0],l[1])
        tnf_dist = getTNFdist(pca_dictionary,l[0],l[1])
        #tnf_dist = getTNFdist(pca_tnf_vectors,c1,c2)
        print(c1 + '\t' + c2 + "\t" + str(dist) + "\t" + str(tnf_dist))



        all_file_names[c1] = True #keep a hash of all filenames
        all_file_names[c2] = True


        '''
        if dist <= threshold:
            if c1 in temp_dictionary.keys():
                temp_dictionary[c2] = temp_dictionary[c1]
            elif c2 in temp_dictionary.keys():
                temp_dictionary[c1] = temp_dictionary[c2]
            else: #neither cols[i] nor cols[j] are in the dictionary yet
                number_of_groups += 1
                if c1 == c2:
                    temp_dictionary[c1] = number_of_groups
                else:
                    temp_dictionary[c1] = number_of_groups
                    temp_dictionary[c2] = number_of_groups


    all_files = all_file_names.keys()
    for filename in all_files:
        #print(filename)
        if filename not in temp_dictionary:
            #print("above file not found")
            unclustered.write(str(filename)+"\n")

    for pair in temp_dictionary.items():
        if pair[1] not in result_dictionary.keys():
            result_dictionary[pair[1]] = []
            print('added ' + str(pair[1]) + ' to result_dictionary')
        result_dictionary[pair[1]].append(pair[0])

    output_metrics.write(str(result_dictionary)+"\n")

    filenames_R1 = {}
    filenames_R2 = {}
    for x in result_dictionary.keys():
        temp_arrayR1 = []
        temp_arrayR2 = []
        for i in result_dictionary[x]:
            if '1.f' in i:
                temp_arrayR1.append(i)
            elif '2.f' in i:
                temp_arrayR2.append(i)

        filenames_R1[x] = ' '.join(sorted(temp_arrayR1))
        filenames_R2[x] = ' '.join(sorted(temp_arrayR2))

    output_metrics.write(str(filenames_R1)+"\n")
    output_metrics.write(str(filenames_R2)+"\n")

    #print('result_dictionary:')
    #print(result_dictionary)

    #FOR INVESTIGATIVE PURPOSES ONLY
    for m in result_dictionary.keys():
        print(Counter(result_dictionary[m]))

    return filenames_R1, filenames_R2
    '''

def main(input_directory,k,s,threshold,cu):

    #get list of all fastq files in directory
    array_of_filenames = glob.glob('./*.fq')

    tnf_vectors = calculateTNF(array_of_filenames)
    pca_tnf_vectors = runPCA(tnf_vectors[0], tnf_vectors[1])
    file_string = ' '.join(array_of_filenames) + ' '

    if cu != 'clusterOnly':
        runMashSketch(file_string,k,s)
        runMashDist(file_string,k,s)
    #concatFilesR1, concatFilesR2 = getClusters('mash_pdist_'+str(k)+'_'+str(s), threshold, pca_tnf_vectors)
    getClusters('mash_pdist_'+str(k)+'_'+str(s), threshold, pca_tnf_vectors)

    '''
    for c in concatFilesR1.keys():
        file_namesR1 = concatFilesR1[c]
        file_namesR2 = concatFilesR2[c]
        cmd1 = 'cat '+file_namesR1 + ' > cluster'+str(c)+'_R1.fastq'
        os.system(cmd1)
        cmd2 = 'cat '+file_namesR2 + ' > cluster'+str(c)+'_R2.fastq'
        os.system(cmd2)
        '''


start = timeit.default_timer()

mash = '~/tools/MASH/mash'

input_directory = sys.argv[1]
iteration = 1
k=12
s=5000
cu = sys.argv[2]
threshold = .95
array_of_filenames=None
os.chdir(input_directory)

pca_tnf_vectors={}

main(input_directory,k,s,threshold,cu)

stop = timeit.default_timer()

print("time to complete: " + str(stop - start))