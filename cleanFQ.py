__author__ = 'KATRINA'
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

