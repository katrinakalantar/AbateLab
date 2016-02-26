__author__ = 'KATRINA'
import sys
import os
import glob
import subprocess


'''
runMashDist.py [directory of fastq files to cluster]
- use this to run MASH sketch and MASH dist for all .fq files in the directory via commands passed to linux

output:
1. .msh sketch files for all corresponding .fq files
3. output file mash_pdist_[k-value]-[s-value] containing the MASH pairwise distances

'''

mash = '~/tools/MASH/mash'
input_directory = sys.argv[1]
iteration = 1
k=12
s=5000
array_of_filenames=None
os.chdir(input_directory)
print(input_directory)

array_of_filenames = glob.glob('./*.f*q')
count = 0

# run MASH sketch (create .msh files from .fq files)
for a in array_of_filenames:
    cmd = mash+' sketch -u -k '+str(k)+' -s '+str(s)+' '+a
    #print(cmd)
    os.system(cmd)
    print('finished runMashSketch')

array_of_sketchFiles = glob.glob('./*.msh')

# run MASH dist (create pairwise distance file from .msh files)
for a in range(len(array_of_sketchFiles)):
    for b in range(a,len(array_of_sketchFiles)):
        print(array_of_sketchFiles[a] + " " + array_of_sketchFiles[b])
        file_string = array_of_sketchFiles[a] + ' ' + array_of_sketchFiles[b]
        if count == 0:
            cmd = mash+' dist -s ' +str(s)+' '+file_string+' > '+'mash_pdist_'+str(k)+'_'+str(s)
            count += 1
        else:
            cmd = mash+' dist -s ' +str(s)+' '+file_string+' >> '+'mash_pdist_'+str(k)+'_'+str(s)
        #print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()