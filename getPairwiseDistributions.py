__author__ = 'KATRINA'

import sys
from matplotlib import pyplot as plt
import numpy as np

#input_file is the mash dist file with pairwise distances
def getPairwiseDistributions(input_file):
    species_dict_1000, species_dict_2500, species_dict_5000, species_dict_7500, species_dict_10000 = {"inter_species":[]},{"inter_species":[]},{"inter_species":[]},{"inter_species":[]},{"inter_species":[]}
    inter_species_dict_1000, inter_species_dict_2500, inter_species_dict_5000, inter_species_dict_7500, inter_species_dict_10000 = {"inter_species":[]},{"inter_species":[]},{"inter_species":[]},{"inter_species":[]},{"inter_species":[]}
    all_species = []

    f = open(input_file,'r').readlines()
    for line in f:
        s = line.strip().split('\t')
        #i1_length = int(s[0].split('/')[1].split('-')[1])
        #i2_length = int(s[1].split('/')[1].split('-')[1])
        #i1_species = s[0].split('/')[1].split('-')[2]
        #i2_species = s[1].split('/')[1].split('-')[2]
        i1_length = int(s[0].split('/')[0].split('-')[1])
        i2_length = int(s[1].split('/')[0].split('-')[1])
        i1_species = s[0].split('/')[0].split('-')[2]
        i2_species = s[1].split('/')[0].split('-')[2]
        i1i2_similarity = s[4].replace('\'','').split('/')
        percent_similarity = float(i1i2_similarity[0])/float(i1i2_similarity[1])

        #this only records if the lengths are both in the same bin
        if i1_length < 1000 and i2_length < 1000:
            if i1_species not in species_dict_1000.keys():
                species_dict_1000[i1_species] = []
                inter_species_dict_1000[i1_species]=[]
            if i1_species == i2_species:
                species_dict_1000[i1_species].append(percent_similarity)
            else:
                inter_species_dict_1000[i1_species].append(percent_similarity)
                species_dict_1000["inter_species"].append(percent_similarity)

        elif i1_length < 2500 and i2_length < 2500:
            if i1_species not in species_dict_2500.keys():
                species_dict_2500[i1_species] = []
                inter_species_dict_2500[i1_species]=[]
            if i1_species == i2_species:
                species_dict_2500[i1_species].append(percent_similarity)
            else:
                inter_species_dict_2500[i1_species].append(percent_similarity)
                species_dict_2500["inter_species"].append(percent_similarity)

        elif i1_length < 5000 and i2_length < 5000:
            if i1_species not in species_dict_5000.keys():
                species_dict_5000[i1_species] = []
                inter_species_dict_5000[i1_species]=[]
            if i1_species == i2_species:
                species_dict_5000[i1_species].append(percent_similarity)
            else:
                inter_species_dict_5000[i1_species].append(percent_similarity)
                species_dict_5000["inter_species"].append(percent_similarity)

        elif i1_length < 7500 and i2_length < 7500:
            if i1_species not in species_dict_7500.keys():
                species_dict_7500[i1_species] = []
                inter_species_dict_7500[i1_species]=[]
            if i1_species == i2_species:
                species_dict_7500[i1_species].append(percent_similarity)
            else:
                inter_species_dict_7500[i1_species].append(percent_similarity)
                species_dict_7500["inter_species"].append(percent_similarity)

        elif i1_length < 10000 and i2_length < 10000:
            if i1_species not in species_dict_10000.keys():
                species_dict_10000[i1_species] = []
                inter_species_dict_10000[i1_species]=[]
            if i1_species == i2_species:
                species_dict_10000[i1_species].append(percent_similarity)
            else:
                inter_species_dict_10000[i1_species].append(percent_similarity)
                species_dict_10000["inter_species"].append(percent_similarity)

    return species_dict_1000, inter_species_dict_1000, species_dict_2500, inter_species_dict_2500, species_dict_5000, inter_species_dict_5000, species_dict_7500, inter_species_dict_7500, species_dict_10000, inter_species_dict_10000




def plotPairwiseDistributions(specD, title):
    species_names = []
    barchart_data = []
    boxplot_data = []
    labels = []
    keysort = sorted(specD.keys())
    for k in keysort:
        species_names.append(k)
        avg = sum(specD[k])/len(specD[k])
        barchart_data.append(avg)
        boxplot_data.append(specD[k])
        labels.append(k)

    figure = plt.figure()
    plt.bar(range(len(barchart_data)), barchart_data, align='center')
    axes = plt.gca()
    axes.set_ylim([0,1])
    plt.xticks(range(len(barchart_data)), labels, rotation=90)# list(barchart_data.keys()), rotation=90)
    plt.title(title)
    plt.tight_layout()
    figure.savefig("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\optimizingForSmallFiles\\100_2\\pairwiseDistibution_"+title+".png",format="png")


    fig2 = plt.figure()
    plt.boxplot(boxplot_data)#, align='center')
    axes = plt.gca()
    axes.set_ylim([0,1])
    plt.xticks(range(len(boxplot_data)), labels, rotation=90) # list(barchart_data.keys()), rotation=90)
    plt.tick_params(axis='x', pad=4)
    plt.title(title)
    plt.tight_layout()
    fig2.savefig("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\optimizingForSmallFiles\\100_2\\pairwiseDistibutionBOX_"+title+".png",format="png")

    plt.close('all')


def comparativePairwiseDistributions(specD, ispecD, title):
    f = open("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\optimizingForSmallFiles\\100_2\\maxmindata.txt",'w')
    for k in specD.keys():
        f.write(k+"\n")
        intra=specD[k]
        inter=ispecD[k]
        if k != 'inter_species':
            f.write('minimum intra species dist: ' + str(min(intra))+"\n")
            f.write('maximum inter species dist: ' + str(max(inter))+"\n")

        bins = np.linspace(0, 1, 200)
        h1 = plt.figure(figsize=(12,5))
        plt.xticks(np.arange(0,1,.05))
        if len(intra)>0:
            plt.hist(intra, bins, alpha=0.5, color="b",label='intra')
        if len(inter)>0:
            plt.hist(inter, bins, alpha=0.5, color="r", label='inter')
        axes = plt.gca()
        axes.set_ylim([0,75])
        plt.legend(loc='upper right')
        plt.title(title+k)
        #plt.show()
        h1.savefig("C:\\cygwin64\\home\\KATRINA\\UCB\\Research\\AbateLab\\optimizingForSmallFiles\\100_2\\pairwiseDistibutionHIST_"+title+k+".png",format="png")

        plt.close()


print('a')
specD1000, ispecD1000, specD2500, ispecD2500, specD5000, ispecD5000, specD7500, ispecD7500, specD10000, ispecD10000 = getPairwiseDistributions(sys.argv[1])

print('begin')

plotPairwiseDistributions(specD1000,'1000')
#plotPairwiseDistributions(specD2500,'2500')
#plotPairwiseDistributions(specD5000,'5000')
#plotPairwiseDistributions(specD7500,'7500')
#plotPairwiseDistributions(specD10000,'10000')

print('second')

comparativePairwiseDistributions(specD1000, ispecD1000, '1000')
print('1000done')
#comparativePairwiseDistributions(specD2500, ispecD2500, '2500')
#print('2500done')
#comparativePairwiseDistributions(specD5000, ispecD5000, '5000')
#comparativePairwiseDistributions(specD7500, ispecD7500, '7500')
#comparativePairwiseDistributions(specD10000, ispecD10000, '10000')


#plotPairwiseDistributions(specD1000,'1000')
