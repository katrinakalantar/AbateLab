__author__ = 'KATRINA'

'''
cleanFQ.py [directory of fastq files to clean]
- use this to remove all reads containing "N" values from .fastq files
- NOTE: this removes reads on a per-file basis, does not consider paired read 
	(ie if only one paired read contains "N" values, paired file output will now contain a read with no pair

	output:
1. fastq files by name [original_fastq_file_name.fastq].o with no reads containing "N" values

'''

import sys
import glob
import os

os.chdir(sys.argv[1])

files = glob.glob('./*.fastq')

for file in files:
    f = open(file,'r').readlines()
    output = open(file+'.o','w')

    count = 0
    l = range(len(f))[0::4]

    for i in l:
        seq = f[i+1]
        if 'N' not in seq:
            output.write(f[i])
            output.write(f[i+1])
            output.write(f[i+2])
            output.write(f[i+3])

