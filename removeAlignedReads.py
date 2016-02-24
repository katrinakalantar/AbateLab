__author__ = 'KATRINA'

'''

removeAlignedReads.py [alignment file (.SAM format from bowtie2)] [root name for paired end .fastq files]

usage: look at the .SAM file and mark all reads which were aligned or unaligned;
remove all aligned reads from the corresponding .fastq files.

- the goal of this is to create .fastq files with 100% novel reads for further analysis and assembly

'''

import sys

alignment_file = sys.argv[1]
fastq_file_root = sys.argv[2]

fq_file1 = open(fastq_file_root+"-1.fq", 'r')    #open paired end fastq files for reading line by line
fq_file2 = open(fastq_file_root+"-2.fq", 'r')
al_file = open(alignment_file, 'r')         #open alignment file for reading line by line

## iterate through alignment file,
## keep an array of all files which were not mapped (to be kept)
reads_to_keep = {}
reads_to_toss = {}
#for line in al_file:
for line in al_file:
    line = line.strip()
    if line[0] == "@":  # Ignore header lines
            continue

    splitLine = line.split("\t")
    name = splitLine[0]             # Sequence name (from .fastq file)
    ref = splitLine[2]              # Reference genome
    seq = splitLine[9]              # Sequence
    flag = int(splitLine[1])        # FLAG
    unmapped = int(bin(flag)[-3])   # Unmapped bit (if 1, read is unmapped)

    if unmapped == 1 and name not in reads_to_toss: #read is not mapped, keep it (only if it doesn't have a mapped pair)
        reads_to_keep[name] = True
    elif unmapped == 0: #this read mapped to a reference, don't keep it in the file
        reads_to_toss[name] = True #NOTE: THIS will prevent keeping a read if one of pairs are unmapped.

#print out read lists for manual checking if necessary
print('reads_to_keep:')
print(str(reads_to_keep))
print('reads_to_toss:')
print(str(reads_to_toss))

#iterate through the R1 .fastq file, create new array with all reads that are being saved
fq_file_output1 = open(fastq_file_root+"-1.fq.m", 'w') #create output fastq file for R1
for line in fq_file1:
    if line[0] == '@':
        seq_name = line.strip()
        if seq_name[1:-2] in reads_to_keep:
            seq = next(fq_file1)
            strand = next(fq_file1)
            qual = next(fq_file1)
            fq_file_output1.writelines(line)
            fq_file_output1.writelines(seq)
            fq_file_output1.writelines(strand)
            fq_file_output1.writelines(qual)


fq_file_output2 = open(fastq_file_root+"-2.fq.m", 'w') #create output fastq file for R2
for line in fq_file2:
    if line[0] == '@':
        seq_name = line.strip()
        if seq_name[1:-2] in reads_to_keep:
            seq = next(fq_file2)
            strand = next(fq_file2)
            qual = next(fq_file2)
            fq_file_output2.writelines(line)
            fq_file_output2.writelines(seq)
            fq_file_output2.writelines(strand)
            fq_file_output2.writelines(qual)


