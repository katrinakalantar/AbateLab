__author__ = 'KATRINA'
import sys
import glob
import os

input_directory = 'C:\\cygwin64\\home\\KATRINA\\UCB\\CourseMaterials\\CS267\\hw1\\benchmarks\\yangsfirstapproach'#sys.argv[1]
os.chdir(input_directory)
array_of_filenames = glob.glob('./*.out')


array_of_all_percentages = []
for file in array_of_filenames:
    print(file)
    f = open(file,'r').readlines()

    all_perc = []
    for line in f:
        if 'Percentage' in line:
            s = line.split(':')
            percentage = s[-1].strip()
            all_perc.append(float(percentage))

    array_of_all_percentages.append(all_perc)
    print(all_perc)





