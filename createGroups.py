__author__ = 'KATRINA'
'''
createGroups.py [directory of fastq files to cluster] [clusterOnly]
- use this to generate clusters for a set of fastq files of the format: name-filelen-1-fastq
- if you already have a mash dist file for the dataset (in the directory specified as arg1)
    type "clusterOnly" in arg2
- else
    type anything else in arg2 (but it does require some jibberish input)

output:
1.clustered files - individual clustered .fq files
2.metrics.o -
3.unclusterd.o - list of files which failed to cluster

*note: any line with the following text: "TODO: switch back for regular usage"
    is intended to deal with fastq files ".fastq" v. ".fq" just uncomment the
    respective lines based on the file format.

'''

import sys
from collections import Counter
import timeit
import os
import glob
import subprocess

#run linux command to create mash sketches for all files
def runMashSketch(file_string,k,s):
    print("inside runMashSketch")
    print(file_string)
    cmd = mash+' sketch -o reference -u -k '+str(k)+' -s '+str(s)+' '+file_string
    print(cmd)
    os.system(cmd)
    print('finished runMashSketch')
    return

#run linux command to create mash distance file
def runMashDist(file_string,k,s):
    print("inside runMashDist")
    cmd = mash+' dist reference.msh -s ' +str(s)+' '+file_string+'> '+'mash_pdist_'+str(k)+'_'+str(s)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)
    print('finished runMashDist')
    return

def getClusters(mash_distanceFile, threshold):
    print("inside getClusters()")
    output_metrics = open("metrics.o",'w')
    unclustered = open("unclustered.o",'w')
    distFile = open(mash_distanceFile,'r')

    result_dictionary = {}; temp_dictionary = {}; all_file_names = {}

    line_count = 0
    number_of_groups = 0

    for line in distFile: #loop through each line in mash distance file
        line_count += 1

        l = line.split('\t')
        c1 = l[0].split('/')[-1]
        c2 = l[1].split('/')[-1]
        dist = 1-(float(l[4].strip().split('/')[0])/float(l[4].strip().split('/')[1]))

        #keep a list of all filenames that have been assigned to a cluster to avoid duplicates
        all_file_names[c1] = True
        all_file_names[c2] = True

        #if the distance < threshold then assign the file to a cluster
        if dist <= threshold:
            print(line)
            #if one file in line already exists in a cluster, add the other filename to that cluster
            if c1 in temp_dictionary.keys():
                temp_dictionary[c2] = temp_dictionary[c1]
            elif c2 in temp_dictionary.keys():
                temp_dictionary[c1] = temp_dictionary[c2]
            else: #neither cols[i] nor cols[j] are in the dictionary yet, begin new cluster
                number_of_groups += 1
                if c1 == c2:
                    temp_dictionary[c1] = number_of_groups
                else:
                    temp_dictionary[c1] = number_of_groups
                    temp_dictionary[c2] = number_of_groups

    print(temp_dictionary)

    #keep output list of files which are not assigned to a cluster
    all_files = all_file_names.keys()
    for filename in all_files:
        if filename not in temp_dictionary:
            unclustered.write(str(filename)+"\n")

    for pair in temp_dictionary.items():
        if pair[1] not in result_dictionary.keys():
            result_dictionary[pair[1]] = []
            print('added ' + str(pair[1]) + ' to result_dictionary')
        result_dictionary[pair[1]].append(pair[0])

    output_metrics.write(str(result_dictionary)+"\n")

    #separate the groups out again based on R1 and R2 paired end files
    filenames_R1 = {}
    filenames_R2 = {}
    for x in result_dictionary.keys():
        temp_arrayR1 = []
        temp_arrayR2 = []
        for i in result_dictionary[x]:
            s = i.split('-')
            #split the name of the file, this depends on how the names are declared so may need to be updated.
            name = s[0]+'-'+s[1] ###TODO: switch bach for regular usage
            #name = s[0]+'-'+s[1]+'-'+s[2]
            temp_arrayR1.append(name+'-1.fastq')  ###TODO: switch bach for regular usage
            temp_arrayR2.append(name+'-2.fastq')
            #temp_arrayR1.append(name+'-1.fq')
            #temp_arrayR2.append(name+'-2.fq')

        '''
        for i in result_dictionary[x]:
            if '1.f' in i:
                temp_arrayR1.append(i)
            elif '2.f' in i:
                temp_arrayR2.append(i)
                '''

        filenames_R1[x] = ' '.join(sorted(temp_arrayR1))
        filenames_R2[x] = ' '.join(sorted(temp_arrayR2))

        print(str(x) + " R1: "+str(filenames_R1))
        print(str(x) + " R2: "+str(filenames_R2))

    output_metrics.write(str(filenames_R1)+"\n")
    output_metrics.write(str(filenames_R2)+"\n")

    #FOR INVESTIGATIVE PURPOSES ONLY
    for m in result_dictionary.keys():
        print(str(m) + str(Counter(result_dictionary[m])))

    return filenames_R1, filenames_R2


#concatonate all paired end read files into one file (prior to running MASH)
def concatonate_PE_reads(array_of_filenames):
    file_names = []
    new_file_array = []
    for i in array_of_filenames:
        s = i.split('-')
        #file_names.append(s[0] + "-" + s[1]+ "-" + s[2]) ###TODO: switch bach for regular usage
        file_names.append(s[0] + "-" + s[1])
    set_names = set(file_names)
    for i in set_names:
        cmd = "cat " + i+"-1.fastq " + i+"-2.fastq > "+i +"-f.fastq"
        #cmd = "cat " + i+"-1.fq " + i+"-2.fq > "+i +"-f.fq"  ###TODO: switch bach for regular usage
        os.system(cmd)
        os.system('wait')
        new_file_array.append(i +"-f.fastq") ###TODO: switch bach for regular usage
        #new_file_array.append(i +"-f.fq")

    return new_file_array



def main(input_directory,k,s,threshold,cu):

    #get list of all fastq files in directory
    #array_of_filenames = glob.glob('./*.fq') ###TODO: switch back for regular usage
    array_of_filenames = glob.glob('./*.fastq')
    #print(array_of_filenames)

    new_file_array = concatonate_PE_reads(array_of_filenames)
    print(new_file_array)

    #convert file array to a string to run through linux command
    file_string = ''
    for f in new_file_array:
        file_string = file_string + f+' '
    file_string = ' '.join(new_file_array) + ' '

    if cu != 'clusterOnly':
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

#declare global variables
start = timeit.default_timer()
mash = '~/tools/MASH/mash'
input_directory = sys.argv[1]
iteration = 1
k=12
s=5000
cu = sys.argv[2]
threshold = .99
array_of_filenames=None
os.chdir(input_directory)

#run main program
main(input_directory,k,s,threshold,cu)

#output total time to console
stop = timeit.default_timer()
print("time to complete: " + str(stop - start))