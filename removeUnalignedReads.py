__author__ = 'KATRINA'

import sys

alignment_file = sys.argv[1]
fastq_file_root = sys.argv[2]

fq_file1 = open(fastq_file_root+"-1.fq", 'r').readlines()    #create paired end fq files - read into arrays
fq_file2 = open(fastq_file_root+"-2.fq", 'r').readlines()
al_file = open(alignment_file, 'r').readlines()         #read alignment file into array

## iterate through alignment file,
## keep an array of all files which were not mapped (to be kept)
reads_to_keep = []
reads_to_toss = []
for line in al_file:
    if line[0] == "@":  # Ignore header lines
            continue

    splitLine = line.split("\t")
    name = splitLine[0]             # Sequence name (from .fastq file)
    ref = splitLine[2]              # Reference genome
    seq = splitLine[9]              # Sequence
    flag = int(splitLine[1])        # FLAG
    unmapped = int(bin(flag)[-3])   # Unmapped bit (if 1, read is unmapped)

    if unmapped == 0 and name not in reads_to_keep and name not in reads_to_toss: #read is not mapped, keep it (only if it isn't already in the list) and doesn't have a mapped pair
        reads_to_keep.append(name)
    elif unmapped == 1:
        reads_to_toss.append(name)  #NOTE: THIS will prevent keeping a read if one of pairs are unmapped.

#print out read lists for manual checking if necessary
print('reads_to_keep:')
print(reads_to_keep)
print('reads_to_toss:')
print(reads_to_toss)

#iterate through the R1 .fastq file, create new array with all reads that are being saved
final_fq_file1 = []
for line_num in range(len(fq_file1)):
    if fq_file1[line_num][0] == '@':
        seq_name = fq_file1[line_num].strip()     # Fastq sequence name
        if seq_name[1:-2] in reads_to_keep:
            final_fq_file1.append(fq_file1[line_num])
            final_fq_file1.append(fq_file1[line_num+1])
            final_fq_file1.append(fq_file1[line_num+2])
            final_fq_file1.append(fq_file1[line_num+3])

#iterate through the R2 .fastq file, create new array with all reads that are being saved
final_fq_file2 = []
for line_num in range(len(fq_file2)):
    if fq_file2[line_num][0] == '@':
        seq_name = fq_file2[line_num].strip()     # Fastq sequence name
        if seq_name[1:-2] in reads_to_keep:
            final_fq_file2.append(fq_file2[line_num])
            final_fq_file2.append(fq_file2[line_num+1])
            final_fq_file2.append(fq_file2[line_num+2])
            final_fq_file2.append(fq_file2[line_num+3])

#write the final output files - .fq.m for "modified"
fq_file_output1 = open(fastq_file_root+"-1.fq.m", 'w') #create output fastq file for R1
fq_file_output1.write(''.join(final_fq_file1))
fq_file_output2 = open(fastq_file_root+"-2.fq.m", 'w') #create output fastq file for R2
fq_file_output2.write(''.join(final_fq_file2))




