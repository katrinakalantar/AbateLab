__author__ = 'KATRINA'
import numpy as np
import pandas as pd
import sys
from collections import Counter
import timeit

'''
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
'''

def main(mash_distanceFile):

    #STRIGHT UP FILE IMPLEMENTATION
    line_count = 0
    threshold = .97
    distFile = open(mash_distanceFile,'r')
    result_dictionary = {}
    temp_dictionary = {}
    metric_dictionary = {}
    number_of_groups = 0
    for line in distFile:
        line_count += 1

        l = line.split('\t')
        c1 = l[0].split('/')[-1]
        c2 = l[1].split('/')[-1]
        dist = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))
        assigned_to_group = False

        if dist <= threshold:
            print(line_count)
            #DICTIONARY IMPLEMENTATION
            if c1 in temp_dictionary.keys():
                temp_dictionary[c2] = temp_dictionary[c1]
            elif c2 in temp_dictionary.keys():
                temp_dictionary[c1] = temp_dictionary[c2]
            else: #neither cols[i] nor cols[j] are in the dictionary yet
                number_of_groups += 1
                print('created a new group: ' + str(number_of_groups))
                if c1 == c2:
                    temp_dictionary[c1] = number_of_groups
                else:
                    temp_dictionary[c1] = number_of_groups
                    temp_dictionary[c2] = number_of_groups

    print(len(temp_dictionary))
    for pair in temp_dictionary.items():
        if pair[1] not in result_dictionary.keys():
            result_dictionary[pair[1]] = []
        result_dictionary[pair[1]].append(pair[0])


    concatonate_dictionary = []
    for x in result_dictionary.keys():
        metric_dictionary[x] = []
        concatonate_dictionary[x] = []
        for i in result_dictionary[x]:
            metric_dictionary[x].append(i.split('-')[-2])
            concatonate_dictionary[x]=result_dictionary[x].join(' ')
    #print(result_dictionary)


            '''
            #LIST IMPLEMENTATION
            for groupk in result_dictionary.keys():
                if c1 in result_dictionary[groupk] and c2 in result_dictionary[groupk]:
                    assigned_to_group = True
                    continue #both already in the result_dictionary, so do nothing for this value
                elif c1 in result_dictionary[groupk] and not c2 in result_dictionary[groupk]:
                    result_dictionary[groupk].append(c2)
                    metric_dictionary[groupk].append(c2.split('-')[-2])
                    assigned_to_group = True
                elif c2 in result_dictionary[groupk] and not c1 in result_dictionary[groupk]:
                    result_dictionary[groupk].append(c1)
                    metric_dictionary[groupk].append(c1.split('-')[-2])
                    assigned_to_group = True
            if assigned_to_group == False: #there was no other group to which this should belong, create new group
                number_of_groups+=1
                print('created a new group: ' + str(number_of_groups))
                result_dictionary[number_of_groups] = []
                metric_dictionary[number_of_groups] = []
                if c1 == c2:
                    result_dictionary[number_of_groups].append(c1)
                    metric_dictionary[number_of_groups].append(c1.split('-')[-2])
                else:
                    result_dictionary[number_of_groups].append(c1)
                    result_dictionary[number_of_groups].append(c2)
                    metric_dictionary[number_of_groups].append(c1.split('-')[-2])
                    metric_dictionary[number_of_groups].append(c2.split('-')[-2])

    '''



    '''


    ##BOOLEAN MATRIX IMPLEMENTATION - PROHIBITIVELY SLOW/LARGE FOR 5000x5000 matrix

    #calculate the distance matrix from mash dist file
    print("creating distance matrix DF")
    df, pairwise_dict = calculateDistanceDF(mash_distanceFile)
    cols = list(df.columns.values)

    print('hello')
    print(df) #NOTE: DF contains similarities
    print(cols)

    result_dictionary = {}
    metric_dictionary = {} #used for investigation of results only
    number_of_groups = 0

    df_matrix = df.as_matrix()
    for i in range(len(df_matrix)):
        for j in range(i,len(df_matrix[i])): #only loop through half of the matrix, since it should be mirrored
            value = df_matrix[i,j] # this is getting each column entry in the row

            assigned_to_group = False


            #LIST IMPLEMENTATION
            if (value): #checks if the value has a 1 in matrix
                for groupk in result_dictionary.keys():
                    if cols[i] in result_dictionary[groupk] and cols[j] in result_dictionary[groupk]:
                        assigned_to_group = True
                        continue #both already in the result_dictionary, so do nothing for this value
                    elif cols[i] in result_dictionary[groupk] and not cols[j] in result_dictionary[groupk]:
                        result_dictionary[groupk].append(cols[j])
                        metric_dictionary[groupk].append(cols[j].split('-')[-2])
                        assigned_to_group = True
                    elif cols[j] in result_dictionary[groupk] and not cols[i] in result_dictionary[groupk]:
                        result_dictionary[groupk].append(cols[i])
                        metric_dictionary[groupk].append(cols[i].split('-')[-2])
                        assigned_to_group = True
                if assigned_to_group == False: #there was no other group to which this should belong, create new group
                    number_of_groups+=1
                    print('created a new group: ' + str(number_of_groups))
                    result_dictionary[number_of_groups] = []
                    metric_dictionary[number_of_groups] = []
                    if cols[i] == cols[j]:
                        result_dictionary[number_of_groups].append(cols[i])
                        metric_dictionary[number_of_groups].append(cols[i].split('-')[-2])
                    else:
                        result_dictionary[number_of_groups].append(cols[i])
                        result_dictionary[number_of_groups].append(cols[j])
                        metric_dictionary[number_of_groups].append(cols[i].split('-')[-2])
                        metric_dictionary[number_of_groups].append(cols[j].split('-')[-2])



            #DICTIONARY IMPLEMENTATION

            if (value):
                if cols[i] in result_dictionary.keys():
                    #if cols[j] not in result_dictionary(): #we can probably assume that it is fine to override this
                    result_dictionary[cols[j]] = result_dictionary[cols[i]]
                elif cols[j] in result_dictionary.keys():
                    result_dictionary[cols[i]] = result_dictionary[cols[j]]
                else: #neither cols[i] nor cols[j] are in the dictionary yet
                    number_of_groups += 1
                    if cols[i] == cols[j]:
                        result_dictionary[cols[i]] = number_of_groups
                    else:
                        result_dictionary[cols[i]] = number_of_groups
                        result_dictionary[cols[j]] = number_of_groups
                '''


    print('result_dictionary:')
    print(result_dictionary)

    #FOR INVESTIGATIVE PURPOSES ONLY
    for m in metric_dictionary.keys():
        print(Counter(metric_dictionary[m]))

start = timeit.default_timer()
main(sys.argv[1])
stop = timeit.default_timer()

print("time to complete: " + str(stop - start))