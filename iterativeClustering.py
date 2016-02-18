import matplotlib
matplotlib.use('Agg')

__author__ = 'KATRINA'

'''
iterativeClustering.py [input directory containing .fastq files to cluster]

usage: this script was used alongside development of clusterData.py to run the hierarchical clustering algorithm
iteratively - recursively calls itself to re-run MASH on clustered files and continue clustering up to 6 iterations.

output:
1. clustered .fastq files in directories corresponding to the iteration
2. linkage_matrix[iteration#].txt - contains the full linkage matrix for each iteration

'''


import sys
import glob
import os
import pandas as pd
from scipy.cluster import hierarchy
from collections import Counter
import seaborn as sns


#input_directory = directory containing .fastq files for each round - will involve copying non-cluster files and grouped files into new directory

def runMashSketch(file_string,k,s,iteration):
    print("inside runMashSketch")
    cmd = mash+' sketch -o reference'+str(iteration)+' -u -k '+str(k)+' -s '+str(s)+' '+file_string
    print(cmd)
    os.system(cmd)
    return

def runMashDist(file_string,k,s,iteration):
    print("inside runMashDist")
    cmd = mash+' dist reference'+str(iteration)+'.msh -s ' +str(s)+' '+file_string+' > '+'mash_pdist_i'+str(iteration)+'_'+str(k)+'_'+str(s)
    print(cmd)
    os.system(cmd)
    return

def calculateDistanceDF(mash_dist_file):
    print("inside calculateDistanceDF")
    pmF = open(mash_dist_file,'r').readlines()
    print(mash_dist_file)
    pairwise_dict = {}
    for line in pmF:
        l = line.split('\t')
        #distances are 1-(fraction of shared kmers in hash)
        pairwise_dict[(l[0].split('/')[-1],l[1].split('/')[-1])] = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))
    ser = pd.Series(list(pairwise_dict.values()),index=pd.MultiIndex.from_tuples(pairwise_dict.keys()))
    df = ser.unstack().fillna(0)
    return [df, pairwise_dict]


def separateClusters(cluster_array, column_headers):
    print("inside sparateClusters")
    max_cluster = max(cluster_array)
    groups = {} #initialize a dictionary of groups
    for i in range(1,max_cluster+1):
        groups[i]=[]
    for x in range(len(cluster_array)): #append each element to dictionary of groups
        groups[cluster_array[x]].append(column_headers[x].split('/')[-1])
    return groups

#input cluster dictionary, plot the composition of each group
def getGroupComposition(cluster_dictionary):
    print("inside getGroupComposition")
    group_id = 0
    iteration_groups_dictionary = {}
    species_counts_dictionary = {}
    for c_list in cluster_dictionary:
        print(str(c_list))
        group_id += 1
        x = cluster_dictionary[c_list]
        species_list = []; file_list = []
        for i in x:
            print(str(i))
            species_list.append(i.split('-')[-2])
            file_list.append(i)
        species_counts = Counter(species_list)
        species_counts_dictionary[group_id] = species_counts
        iteration_groups_dictionary[group_id] = file_list
    return iteration_groups_dictionary, species_counts_dictionary #species_groups_dictionary

def combineGroupFiles(array_of_filenames, species, i, iteration):
    files_to_combine = ' '.join(array_of_filenames)
    cmd = 'cat '+files_to_combine+' > iter'+str(iteration)+'_group'+str(i)+'-'+species+'-c.fq'
    os.system(cmd)
    return

mash = '~/tools/MASH/mash'

clusters = None
uncategorized_clusters_exist = True
input_directory = sys.argv[1]
iteration = 1
k=16
s=400
array_of_filenames=None
os.chdir(input_directory)
valid_merge=0
invalid_merge=0
false_negative=0
false_positive=0

count1 = 0


def main(input_directory, array_of_filenames, uncategorized_clusters_exist,iteration):

    print("inside main")
    print("ITERATION #"+str(iteration))

    if iteration > 6:
        print("running time has exceeded three recursive iterations")
        return

    #get list of all fastq files in directory
    if len(array_of_filenames) < 1:# == None:
        array_of_filenames = glob.glob('./*.fq')

    #convert file array to a string to run through linux command
    file_string = ''
    for f in array_of_filenames:
        file_string = file_string + f+' '

    if uncategorized_clusters_exist or count1 < 2:
        print("inside loop - uncategorized_clusters_exist")
        uncategorized_clusters_exist = False #reset boolean to False, will set to True again if we find more uncategorized clusters
        runMashSketch(file_string,k,s,iteration)
        runMashDist(file_string,k,s,iteration)

        #calculate the distance matrix from mash dist file
        print("creating distance matrix DF")
        df, pairwise_dict = calculateDistanceDF('mash_pdist_i'+str(iteration)+'_'+str(k)+'_'+str(s))
        cols = list(df.columns.values)

        #run clustering on the distance matrix
        print("begin running clustering on distance matrix")
        linkage_matrix = hierarchy.linkage(df, method='ward', metric='euclidean')
        f=open('linkage_matrix'+str(iteration)+'.txt','w')
        f.write(str(linkage_matrix))

        ###adjusting maxclust based on iteration
        m = 75-5*iteration
        print("maxclust: "+str(m))

        clusters = hierarchy.fcluster(linkage_matrix,m,'maxclust') ##THIS IS THE LINE I WANT TO USE WITH KNOWN CLUSTER NUMBERS
        #clusters=hierarchy.fcluster(linkage_matrix,1.5+(.3*iteration),'inconsistent',depth=9)
        print(clusters)

        cluster_dictionary = separateClusters(clusters, cols)
        groups, species_counts = getGroupComposition(cluster_dictionary)

        ##1/20 getting within_cluster distance:
        avg_cluster_distances = []

        for clust in groups:
            dist_array = []
            for a in groups[clust]:
                for b in groups[clust]:
                    if a != b:
                        d = pairwise_dict[(a,b)]
                        dist_array.append(d)
            try:
                avg = sum(dist_array)/len(dist_array)   #SOME CLUSTERS HAVE ONLY 1 ENTRY IN THEM -- NOT SURE WHY THIS IS HAPPENING, BUT IT CAUSES 0 length dist[] because a==b
            except:
                avg = 0
            avg_cluster_distances.append(avg)
        #print(avg_cluster_distances)

        #sns_plot3 = sns.clustermap(df, row_linkage = linkage_matrix, col_linkage = linkage_matrix, method='ward', cbar = True, cmap = "summer_r", figsize=(8,8))
        #fig3 = sns_plot3.savefig("sns_clustermap_output_"+str(iteration)+".png")


        next_file_names = []
        for i in range(len(species_counts)):
            #if len(species_counts[i+1]) == 1: #this is a pure species file, combine all groups  ##1/22 changing this to distance metric
            if avg_cluster_distances[i] < 0.99:   #CAN"T TELL IF THIS SHOULD BE i or i+1 - actually, should be i because its an array and the others are dictionaries with int keys
                print("~valid cluster by avg dist~")
                #print("cluster contains single species")
                global valid_merge
                global false_positive
                if len(species_counts[i+1]) == 1:
                    valid_merge+=1
                else:
                    false_positive +=1
                combineGroupFiles(groups[i+1],list(species_counts[i+1].keys())[0],i+1,iteration)
                next_file_names.append('iter'+str(iteration)+'_group'+str(i+1)+'-'+list(species_counts[i+1].keys())[0]+'-c.fq')
            else: #this file contains many different species, do not combine + add all file names to list
                global false_negative
                global invalid_merge
                if len(species_counts[i+1]) == 1:
                    false_negative+=1
                else:
                    invalid_merge+= 1

                for g in groups[i+1]:
                    next_file_names.append(g)
                uncategorized_clusters_exist = True
                if count1 == 0:
                    count = 1
        print(next_file_names)

        iteration += 1
        main(input_directory,next_file_names,uncategorized_clusters_exist,iteration)


print("hello world")
main(input_directory, {}, uncategorized_clusters_exist,iteration)

print("invalid merge count: " + str(invalid_merge))
print("valid merge count: " + str(valid_merge))
print("false positives:" + str(false_positive))
print("false negatives:" + str(false_negative))
