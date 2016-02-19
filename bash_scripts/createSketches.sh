#!/usr/bin

#usage: bash globalMash_tophits.sh test1
#note - no trailing slash after the root directory name

echo !!BEGINNING SKETCHING!!
date

mash=~/tools/MASH/mash
root_directory=$1
input_directory="$1/*"
#refseq_hash=~/tools/MASH/RefSeqSketches.msh

mkdir "$root_directory/sketches"
files="$root_directory/*.fq"

for f in $files
do
echo $f
mash sketch -u -k 16 -s 400 $f
done

cd $root_directory
mv *.msh ./sketches/
cd -

echo !!SKETCHING COMPLETED!!
date


