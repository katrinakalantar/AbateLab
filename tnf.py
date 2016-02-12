import sys
import itertools
import glob
import os
from matplotlib.mlab import PCA
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from matplotlib import cm
from sklearn.lda import LDA


def getKmerFreqs(input_file):
    f = open(input_file,'r').readlines()

    alphabet = 'ACGT'
    kmers = [''.join(i) for i in itertools.product(alphabet, repeat = 4)]

    kmerD = {} #initialize kmer dictionary
    for k in kmers:
        kmerD[k] = 0

    total_kmer_count = 0
    for line_num in range(len(f)):
        if line_num%4 == 1:
            line = f[line_num]
            for i in range(len(line)-4):
                total_kmer_count += 1
                kmerD[line[i:i+4]] += 1

    for k in kmerD:
        kmerD[k] = kmerD[k]/total_kmer_count

    s = sorted(kmerD)
    kmer_result_vector = []
    for value in s:
        kmer_result_vector.append(kmerD[value])

    return kmer_result_vector


start = time.time()
input_dir = sys.argv[1]
os.chdir(input_dir)
array_of_filenames = glob.glob('./*.fq')

all_kmer_vectors = {}
all_kmer_vectors_array = []
labels = []
for input_file in array_of_filenames:
    labels.append(str(input_file).split('-')[2])
    v = getKmerFreqs(input_file)
    #print(v)
    all_kmer_vectors_array.append(np.asarray(v))
    all_kmer_vectors[input_file] = v

#print(labels)

print('finished kmer frequency calculations: ' + str(time.time()-start))
#print(all_kmer_vectors)


def runPCA(all_kmer_vectors_array):

    results = PCA(np.asarray(all_kmer_vectors_array))

    print(results)
    print(results.fracs)
    print(results.Y)

    print(len(results.fracs))
    print(sum(results.fracs))
    print(len(results.Y))

    x=[]
    y=[]
    z=[]
    c=[]
    for item in results.Y:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])
        c.append(item[3])

    fig1 = plt.figure()
    ax = Axes3D(fig1)
    pltData = [x,y,z,c]
    ax.scatter(pltData[0],pltData[1],pltData[2], c=pltData[3])
    for i in range(len(x)):
        if i %20 == 0: #only write every 1/20 labels
            ax.text(x[i],y[i],z[i],  '%s' % (str(labels[i])), size=5, zorder=1, color='k')
        print(str(c[i]) + str(labels[i]))

    # make simple, bare axis lines through space:
    xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
    yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
    zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.

    # label the axes
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title("PCA analysis of sim2_1000fq_0129_ocean TNF data")
    plt.show() # show the plot

    return

print('before runPCA')
runPCA(all_kmer_vectors_array)
print('after runPCA')

def runLDA(all_kmer_vectors_array,labels):
    sklearn_lda = LDA(n_components=4)
    X = np.array(all_kmer_vectors_array)
    y = np.array(labels)
    X_lda_sklearn = sklearn_lda.fit_transform(X,y)
    print(X_lda_sklearn)
    return X_lda_sklearn


def plot_scikit_lda(LDAresult):

    x=[]
    y=[]
    z=[]
    c=[]
    for item in LDAresult:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])
        c.append(item[3])

    fig1 = plt.figure()
    ax = Axes3D(fig1)
    pltData = [x,y,z,c]
    ax.scatter(pltData[0],pltData[1],pltData[2], c=pltData[3] )
    for i in range(len(x)):
        if i %20 == 0: #only write every 1/20 labels
            ax.text(x[i],y[i],z[i],  '%s' % (str(labels[i])), size=5, zorder=1, color='k')
        #print("["+str(x[i])+","+str(y[i])+","+str(z[i])+","+str(c[i])+"]"+ str(labels[i]))

    '''
    for label,x,y,z in zip(labels,x,y,z):
        ax.annotate(
            label, xyz=(x,y,z)
        )
        '''

    # make simple, bare axis lines through space:
    xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0), (0,0)) # 2 points make the x-axis line at the data extrema along x-axis
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r') # make a red line for the x-axis.
    yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])), (0,0)) # 2 points make the y-axis line at the data extrema along y-axis
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r') # make a red line for the y-axis.
    zAxisLine = ((0, 0), (0,0), (min(pltData[2]), max(pltData[2]))) # 2 points make the z-axis line at the data extrema along z-axis
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r') # make a red line for the z-axis.

    # label the axes
    ax.set_xlabel("LD1")
    ax.set_ylabel("LD2")
    ax.set_zlabel("LD3")
    ax.set_title("LDA analysis of sim2_1000fq_0129_ocean TNF data")
    plt.show() # show the plot





    '''
    fig1 = plt.figure()
    ax = Axes3D(fig1)
    for label,marker,color in zip(
        range(1,21),('^','^','^', 's','s','s', 'o','o','o','x','x','x','8','8','8','1','1','1','+','+'),('blue', 'red', 'green','blue', 'red', 'green','blue', 'red', 'green','blue', 'red', 'green','blue', 'red', 'green','blue', 'red', 'green','blue', 'red')):

        plt.scatter(x=X[:,0][labels == label]*mirror,
                y=X[:,1],z=X[:,2][labels == label],
                marker=marker,
                color=color,
                alpha=0.5,
                label=label_dict[label]
                )

    ax.set_xlabel("LD1")
    ax.set_ylabel("LD2")
    ax.set_zlabel("LD3")

    leg = plt.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    plt.title(title)

    # hide axis ticks
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
            labelbottom="on", left="off", right="off", labelleft="on")

    # remove axis spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.grid()
    #plt.tight_layout
    '''
    plt.show()


label_dict = {1: 'Setosa', 2: 'Versicolor', 3:'Virginica',4:'Virginica',5:'Virginica',6:'Virginica',7:'Virginica',8:'Virginica',9:'Virginica',10:'Virginica',11:'Virginica',12:'Virginica',13:'Virginica',14:'Virginica',15:'Virginica',16:'Virginica',17:'Virginica',18:'Virginica',19:'Virginica',20:'Virginica'}
LDAresult = runLDA(all_kmer_vectors_array,labels)
plot_scikit_lda(LDAresult)