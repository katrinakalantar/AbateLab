__author__ = 'KATRINA'

'''

krakenDistributions.py [input directory] [kraken file]
- take in the kraken file (obtained from running ~/scripts/fullKraken.sh > filename) (on AbateLab server)
file contents should look like:
"
1000-143-DesulfotaleaPsychrophilaLSv54-1.fq

Processed 143 sequences (21450 bp) ...
143 sequences (0.02 Mbp) processed in 0.012s (695.7 Kseq/m, 104.35 Mbp/m).
  140 sequences classified (97.90%)
  3 sequences unclassified (2.10%) "

- assumes that there are *.fq-kr-ov files in the input directory
- will use the input file kraken file alongside the *.fq-kr-ov to determine the dominant purity

output:
kraken_purity.txt - output file containing filename and purity information

*note: the plotting functions should work, but we don't need a million pie charts,
 so they are commented out by default

'''

import sys
import glob
import os
import matplotlib.pyplot as plt
import operator

#this function is not in use currently
def createPieChart(labels, values):
    plt.pie(values, labels=labels)
    plt.show()

#open each .fq-kr-ov file and load the distribution of classified reads and counts
#append this distribution to the summary dictionary - used to generate the .tsv purity file later
def loadDistributions(array_of_filenames, summary_dictionary):
    massive_dictionary = {}
    total_of_all_reads = {}
    for file in array_of_filenames:
        print(file)
        f = open(file,'r').readlines()
        per_file_total = 0
        species_dict = {}
        for line in f:
            sp = line.split()
            if len(sp)<2: # this is the line corresponding to unclassified
                species_dict['low_classification'] = int(sp[0])
            else:
                species_dict[sp[1]] = int(sp[0])

        print(file.split('/')[-1].split('.')[0])
        summary_dictionary[file.split('/')[-1].split('.')[0]].append(species_dict)
        sorted_species_dict = sorted(species_dict)

        labels = []; values = []
        for v in sorted_species_dict:
            labels.append(v)
            values.append(species_dict[v])
            if v in total_of_all_reads.keys():
                total_of_all_reads[v] += species_dict[v]
            else:
                total_of_all_reads[v] = species_dict[v]

        massive_dictionary[file.split('\\')[-1].split('.')[0]] = [labels, values]
    return massive_dictionary, total_of_all_reads, summary_dictionary


#write the purity data obtained in logic above into a file
def writeOutputFile(summary_dictionary):
    f = open("purity_output.txt",'w')
    for i in summary_dictionary:
        output_string = [i]
        statistics = summary_dictionary[i]
        classified = statistics[0]
        unclassified = statistics[1]

        if (len(statistics) > 2):
            other = statistics[2]
            if classified > unclassified:
                max_species_value = max(other.values())
                total_values = sum(other.values())
                x = max_species_value*classified / total_values
                for k in other:
                    if other[k] == max_species_value:
                        max_species = k
                output_string.append(str(x)) ##TODO !! THIS SHOULD BE A NORMALIZED DECIMAL VALUE IDEALLY
                output_string.append(max_species)

            else:
                output_string.append(str(unclassified))
                output_string.append('unclassified')

            o_sorted = reversed(sorted(other.items(), key=operator.itemgetter(1)))
            for i in o_sorted:
                v = i[0]+":"+str(i[1])
                output_string.append(v)
            f.write('\t'.join(output_string)+"\n")

        else:
            if(unclassified > classified):
                output_string.append(str(unclassified))
                output_string.append('unclassified')
                f.write('\t'.join(output_string)+"\n")


#MAIN SCRIPT STARTS HERE
input_dir = sys.argv[1]
os.chdir(input_dir)
os.system('pwd')
array_of_filenames = glob.glob('./*-ov')

out_file_name=sys.argv[2]
out_file = os.path.join(input_dir,out_file_name)
f = open(out_file_name,'r').readlines()

summary_dictionary = {}
percentage_classified = []
percentage_unclassified = []
bargroups = []
for line in f:
    summary_line = []
    if '.fastq' in line: #note if file names are *.fq this will need to be updated
        bargroups.append(line.strip().split('.')[0])
    if 'sequences classified' in line:
        percent_c = float(line.split('(')[1].split('%')[0])
        percent_u = 100-percent_c
        percentage_classified.append(percent_c)
        percentage_unclassified.append(percent_u)

#add bargroup name, percent classified, and percent unclassified to the summary dictionary
for i in range(len(bargroups)):
    summary_dictionary[bargroups[i]] = [percentage_classified[i],percentage_unclassified[i]]

#create a histogram of percent unclassified reads
'''
fig1 = plt.figure()
plt.hist(percentage_unclassified,bins=40,color='#00bfff')
plt.title('distribution of the percent of unclassified reads')
plt.show()
'''

distributions = loadDistributions(array_of_filenames, summary_dictionary)
stats = distributions[2] ##this will be the summary array which I will want to post-process into .tsv
md = distributions[0] #for plotting the pie charts below

'''#create pie charts for all barcode groups
for x in md:
    print(x)
    print(md[x][0])
    print(md[x][1])
    createPieChart(md[x][0],md[x][1])
'''
totals = distributions[1]
#createPieChart(totals.keys(),list(totals.values()))

writeOutputFile(summary_dictionary)