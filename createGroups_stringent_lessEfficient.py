__author__ = 'KATRINA'
import numpy as np
import pandas as pd
import sys
from collections import Counter
import timeit
import os
import glob



def runMashSketch(file_string,k,s):
    print("inside runMashSketch")
    cmd = mash+' sketch -o reference -u -k '+str(k)+' -s '+str(s)+' '+file_string
    print(cmd)
    os.system(cmd)
    return

def runMashDist(file_string,k,s):
    print("inside runMashDist")
    cmd = mash+' dist reference.msh -s ' +str(s)+' '+file_string+' > '+'mash_pdist_'+str(k)+'_'+str(s)
    print(cmd)
    os.system(cmd)
    return


def getClusters(mash_distanceFile, threshold):
    print("inside getClusters()")
    output_metrics = open("metrics.o",'w')
    distFile = open(mash_distanceFile,'r')

    result_dictionary = {}
    temp_dictionary = {}
    #metric_dictionary = {}

    line_count = 0
    number_of_groups = 0

    for line in distFile:
        line_count += 1
        assigned_to_group = False

        l = line.split('\t')
        c1 = l[0].split('/')[-1]
        c2 = l[1].split('/')[-1]
        dist = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))

        if dist <= threshold:
            if c1 in temp_dictionary.keys():
                if '928-808-ProchlorococcusMarinusCCMP1375' in c1:
                    print('-----')
                    print(c1)
                    print(c2)
                    print(c1 in temp_dictionary.keys())
                    print(c2 in temp_dictionary.keys())
                    print(temp_dictionary[c1])
                if c2 not in temp_dictionary.keys():
                    temp_dictionary[c2] = temp_dictionary[c1]
            elif c2 in temp_dictionary.keys():
                if '928-808-ProchlorococcusMarinusCCMP1375' in c2:
                    print('-----')
                    print(c1)
                    print(c2)
                    print(c1 in temp_dictionary.keys())
                    print(c2 in temp_dictionary.keys())
                    print(temp_dictionary[c2])
                if c1 not in temp_dictionary.keys():
                    temp_dictionary[c1] = temp_dictionary[c2]
            else: #neither cols[i] nor cols[j] are in the dictionary yet
                number_of_groups += 1
                print('created a new group: ' + str(number_of_groups))
                output_metrics.write('created a new group: ' + str(number_of_groups)+"\n")
                output_metrics.write(str(temp_dictionary)+"\n")
                print('c1: ' + str(c1))
                print('c2: ' + str(c2))
                print('hello')
                if c1 == c2:
                    print('A')
                    temp_dictionary[c1] = number_of_groups
                    print(temp_dictionary[c1])
                else:
                    print('B')
                    temp_dictionary[c1] = number_of_groups
                    temp_dictionary[c2] = number_of_groups
                    print(temp_dictionary[c1])
                    print(temp_dictionary[c2])


    #print(str(temp_dictionary))

    #print(len(temp_dictionary))
    for pair in temp_dictionary.items():
        if pair[1] not in result_dictionary.keys():
            result_dictionary[pair[1]] = []
            print('added ' + str(pair[1]) + ' to result_dictionary')
        result_dictionary[pair[1]].append(pair[0])

    output_metrics.write(str(result_dictionary))

    concatonate_dictionary_R1 = {}
    concatonate_dictionary_R2 = {}
    filenames_R1 = {}
    filenames_R2 = {}
    for x in result_dictionary.keys():
        #metric_dictionary[x] = []
        concatonate_dictionary_R1[x] = []
        concatonate_dictionary_R2[x] = []
        for i in result_dictionary[x]:
            if '1.f' in i:
                concatonate_dictionary_R1[x].append(i)
            elif '2.f' in i:
                concatonate_dictionary_R2[x].append(i)

        filenames_R1[x] = ' '.join(concatonate_dictionary_R1[x])
        filenames_R2[x] = ' '.join(concatonate_dictionary_R2[x])

    output_metrics.write(str(filenames_R1))
    output_metrics.write(str(filenames_R2))

    #print('result_dictionary:')
    #print(result_dictionary)

    #FOR INVESTIGATIVE PURPOSES ONLY
    #for m in metric_dictionary.keys():
    #    print(Counter(metric_dictionary[m]))

    return filenames_R1, filenames_R2


def main(input_directory,k,s,threshold):

    #get list of all fastq files in directory
    array_of_filenames = glob.glob('./*.fq')

    #convert file array to a string to run through linux command
    file_string = ''
    for f in array_of_filenames:
        file_string = file_string + f+' '


    runMashSketch(file_string,k,s)
    runMashDist(file_string,k,s)
    concatFilesR1, concatFilesR2 = getClusters('mash_pdist_'+str(k)+'_'+str(s), threshold)

    for c in concatFilesR1.keys():
        file_namesR1 = concatFilesR1[c]
        file_namesR2 = concatFilesR2[c]
        cmd1 = 'cat '+file_namesR1 + ' > cluster'+str(c)+'_R1.fastq'
        os.system(cmd1)
        cmd2 = 'cat '+file_namesR2 + ' > cluster'+str(c)+'_R2.fastq'
        os.system(cmd2)


start = timeit.default_timer()

mash = '~/tools/MASH/mash'

input_directory = sys.argv[1]
iteration = 1
k=16
s=400
threshold = .97
array_of_filenames=None
os.chdir(input_directory)


main(input_directory,k,s,threshold)

stop = timeit.default_timer()

print("time to complete: " + str(stop - start))