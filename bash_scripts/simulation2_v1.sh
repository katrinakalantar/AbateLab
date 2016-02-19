#/bin/bash
#usage: bash simulation1_f.sh simulation_name 
#uses ART to simulates illumina reads for random genomes, with randomized lengths

echo !!BEGINNING SIMULATION!!
date

art=~/simulation/art_bin_ChocolateCherryCake/art_illumina
reference_directory=~/simulation/ocean_reference_files/
simulation_name=$1
mkdir $simulation_name

for i in {1..1000}
do
echo $i

##GET RANDOM INT for file length
len=$(python -S -c "import random; print random.randrange(100,102)")

##CHOOSE A REFERENCE FILE at random from reference directory
cd $reference_directory
reference_file="$(find . -type f | shuf -n 1)"
reference_file2=${reference_file:2}
ref_name="${reference_file2%%.*}"
cd ../scripts
echo $ref_name

#RECORD reference file used fro calculating distribution of genomes
echo $ref_name >> reference_genomes_used.txt
a="$ref_name-"
b="$i-"
c="$len-"
output_file_name="$b$c$a"

#SIMULATION of PE reads; 150bp; mean fragment sizes 500; standard deviation 10
#echo $len is len parameter
$art -i "$reference_directory$reference_file" -o "$output_file_name" -ss MS -l 150 -c $len -p -m 500 -s 10 -q 
#wait

#ATTEMPTING TO GET THE ACTUAL NUMBER OF READS FROM THE -c PARAM - SOMETHING GNG WRONG HERE
#the debug statements are only for investigation, remaining lines are required for cut in final step
ext1=1.fq
ext2=2.fq
output1=$output_file_name$ext1
output2=$output_file_name$ext2
echo here1

actual_len="$(wc -l $output1|awk '{print $1}')"
count=4

actual_reads=$((actual_len/count))
echo here2
echo $((actual_len/count)) >> distribution_of_actual_reads_per_fasta.txt

#RECORD number of reads EXPECTED in each file
echo $len >> distribution_of_reads_per_fasta.txt 

#IF THE actual length of the file doesn't match expected
#due to something weird going on with -c parameter
#then cut the output to only the originally intended length
if [ $actual_reads != $len ]; then
	echo here3
	echo WARNING actual reads do not match expected
	new_len=$((len*4))
	head -$new_len $output1 > temp
	cat temp > $output1
	head -$new_len $output2 > temp
	cat temp > $output2
fi


#MOVE files to final directory
mv *.aln $simulation_name
mv *.fq $simulation_name

done

mv *.txt $simulation_name

echo !!COMPLETED SIMULATION!!
date
