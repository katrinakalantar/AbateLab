# AbateLab
Gapped Genomes Clustering 
work done in the Abate Lab at UCSF (Winter 2016 Rotation)

# Most Important Scripts:
createGroups.py
krakenDistributions.py



### cleanFQ.py
cleanFQ.py [directory containing .fastq files to clean]

usage: remove all reads containing "N" values from .fastq files
- NOTE: this removes reads on a per-file basis, does not consider paired read 
	(ie if only one paired read contains "N" values, paired file output will now contain a read with no pair
output:
1. fastq files by name [original_fastq_file_name.fastq].o with no reads containing "N" values

### clusterByBoolMatrix.py
clusterByBoolMatrix.py [mash distance file]

usage: OUTDATED - was originally used to cluster via pairwise distances before implementing threshold-based clustering approach

### clusterData.py
clusterData.py [input data matrix] [input refseq taxons] [output file name]

usage: this script was used to cluster data based on MASH distance to the RefSeq database. 
Creates a dendrogram from hierarchical clustering (via python fastcluster module) of the refseq MASH data.

- NOTE: prohibitively slow for large datasets; this approach was abandoned after pursuing pairwise method
output:
1. [output file name]_dendroMatrix.txt - saves the input to the dendrogram after running fastcluster
2. [output file name]_dendro.png - saves the .png image of the dendrogram


### clusterData2.py
clusterData2.py [pairwise MASH distance file] [output file name]

usage: this does a full clustering (without actually concatonating files) based on the hierarchical clustering
and linkage graph. This script will fail for input sizes > 2000 fastq files - linkage function and hierarchical
clustering are computationally intensive.

- NOTE: this approach was abandoned due to difficulty setting a cutoff on the dendrogram for clustering

output:
1. outputs several files and images for analysing the species composition of each cluster 
2. [output file name]sns_clustermap_output.png - heatmap of the clustered data


### createGroups.py
createGroups.py [directory of fastq files to cluster] [clusterOnly]
- use this to generate clusters for a set of fastq files of the format: name-filelen-1-fastq
- if you already have a mash dist file for the dataset (in the directory specified as arg1)
    type "clusterOnly" in arg2
- else
    type anything else in arg2 (but it does require some jibberish input)

output:
1.clustered files - individual clustered .fq files
2.metrics.o -
3.unclusterd.o - list of files which failed to cluster

*note: any line with the following text: "TODO: switch back for regular usage"
    is intended to deal with fastq files ".fastq" v. ".fq" just uncomment the
    respective lines based on the file format.

### createGroups_stringent_lessEffective
OUTDATED - precursor file to createGroups.py (above) used to test parameters for k, s, and threshold values

### getPairwiseDistributions.py
getPairwiseDistributions.py [pairwise MASH dist file]

usage: this script was used to plot histograms of inter- and intra- species similarity and set 
threshold values for differing file sizes (reads/fastq)
*entirely for investigative purposes

output: 
per-species, per-file-size-category plots of inter- and intra-species distance (outputs LOTS of plots)
	
	
### investigatingVectors
OUTDATED - this was used to understand the composition of vectors resulting from MASH distance comparison to RefSeq DB
this methods was abandoned in pursuit of the pairwise MASH distance approach
	
	
### iterativeClustering.py
iterativeClustering.py [input directory containing .fastq files to cluster]

usage: this script was used alongside development of clusterData.py to run the hierarchical clustering algorithm
iteratively - recursively calls itself to re-run MASH on clustered files and continue clustering up to 6 iterations.

output:
1. clustered .fastq files in directories corresponding to the iteration
2. linkage_matrix[iterationNumber].txt - contains the full linkage matrix for each iteration

	
### krakenDistributions.py
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


### loadData.py
loadData.py [input data matrix] [input refseq taxons] [output file name]

OUTDATED - January 8th 2015
This script does the following:
1. Imports a file of mash-refseq hit vectors
2. Imports a file of refseq taxonomic classifications
3. Removes all rows containing zero values for all fastq files
4. Removes the corresponding refseq taxonomic classifications
5. Writes the new data matrix and refseq data to files


### pairwiseAnalysis.py
pairwiseAnalysis.py [MASH distance file] [output file name]

usage: create a cluster map (heatmap) of the data in the MASH pairwise distance file; used this to visualize the potential for a sample to cluster

output:
1. [output file name]sns_clustermap_output.png


### viewAsGraph.py
viewAsGraph.py [MASH distance file]

usage: use this to view the distance file as a graph - 
- NOTE: cannot run on abatelab computer bc of package dependencies; only works for small sets of files when run locally

output:
1. graph.png - graph representation of the MASH distance file

