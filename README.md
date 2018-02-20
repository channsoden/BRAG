# BRAG [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/channsoden/BRAG/master)

Break Rates Across the Genome (BRAG) is a program to estimate the rearrangement break rate at all sites of a genome, given a tree of relatives and pairwise alignments between the genome of interest (the reference) and all other leaves in the tree.

## Usage

### Mandatory Arguments
--tree  
Newick formatted phylogenetic tree should have the reference and all other genomes to be compared, and only the other genomes to be compared, including one outgroup.

--reference  
The name of the reference, as labelled in the tree.

--outgroup  
The name of the outgroup, as labelled in the tree.

--segment_files  
List (separated by spaces) of OSFinder-formatted pairwise alignments between the reference and all queries in the tree.

--seqs_files  
List (separated by spaces) of files containing the paths to the genomes used in each pairwise alignment (segment_files), as produced by Murasaki. These files should be in the same order as the alignments. Each should have two lines. The first line should be the path to the fasta-formatted genome that appears first in the corresponding alignment (lefthand columns). The second line should be the path to the fasta-formatted genome that appears second in the corresponding alignment (righthand columns). The seqs file names should match the following pattern:

PREFIX.SHORTNAME1.SHORTNAME2.*

Where PREFIX, SHORTNAME1, and SHORTNAME2 do not contain any "." characters, and SHORTNAME1 AND SHORTNAME2 are identifiers for the reference and a query, one of which should match the --reference option. The * can represent any or no string, including "." characters.

### Optional Arguments
--output  
A name to be used as a prefix for output files. Default "rearrangement_analysis".

--centromeres  
A tab separated table of centromere coordinates. There should be 3 columns.
1. The index of the scaffold, starting with 0, in the order of the sequences
   in the reference genome.
2. The start coordinate of the centromere on that scaffold.
3. The stop coordinate of the centromere on that scaffold.  
Comments after "#" will be ignored.

--tracks  
A tab separated table of extra data tracks to plot along the reference genome, and with which to perform simple linear regression against the break rate.The first line should be a header containing column labels. There should be 2 columns labeled 'start' and 'end', which contain the start and end coordinates of windows along the reference genome using a half-closed interval: [start, end). The window
size and step of these windows should match the window size and step arguments given to BRAG.

E.G. with the default window size of 35kbp and a 7kbp step:

start	end	GC  
0	35000	0.4301714286  
7000	42000	0.4380285714  
14000	49000	0.4396571429  
21000	56000	0.4472857143  
Missing data can be left blank.

--step  
The step size for sliding window analysis (AKA rectangular smoothing).

--window_size  
The window size for sliding window analysis (AKA rectangular smoothing).

--nthreads  
The number of threads to use for some steps. Default 1.

## Sample Analysis

To run the sample analysis, from the sample_data directory, simply run:

```bash
bash sample_analysis.sh
```

This analysis can be run on MyBinder.org from your browser at the following link:

https://mybinder.org/v2/gh/channsoden/BRAG/master


