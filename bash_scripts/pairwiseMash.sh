#!/bin/bash
#TRY TO USE THE PAIRWISE COMPARISON METHOD
#http://mash.readthedocs.org/en/latest/tutorials.html#pairwise-comparisons-with-compound-sketch-files

#default uses kmer size = ; sketch size = 1000 

echo !!STARTING PAIRWISE MASH COMPARISON SCRIPT!!
date

root_directory=$1 #name of directory containing fastq files (omit trailing slash)
mash=~/tools/MASH/mash
files="$root_directory/*.fq"
#k = $2 ; s = $3


echo !!GENERATING pairwise reference file!!
date
if [[$# -eq 3]]; then
$mash sketch -o reference -u -k $2 -s $3 $files #fq1 fq2 fq3...fqN
du ~/sandbox/test1_pairwise/reference.msh
else #run with default parameters
$mash sketch -o reference -u $files #fq1 fq2 fq3...fqN
fi

echo !!CALCULATING pairwise distances!!
date
if [[$# -eq 3]]; then
$mash dist reference.msh -s $3 $files > "$root_directory/mash_pairwiseDistances-$2-$3" #fq1 fq2 fq3...fqN # THIS contains ALL pairwise distances in out output file
du "$root_directory/mash_pairwiseDistances-$2-$3"
else #run with default parameters
$mash dist reference.msh $files > "$root_directory/mash_pairwiseDistances"
fi
mv reference.msh $root_directory

echo !!COMPLETED pairwise distance script!!
date
