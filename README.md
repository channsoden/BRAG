# BRAG

Break Rates Across the Genome (BRAG) is a program to estimate the rearrangement break rate at all sites of a genome, given a tree of relatives and pairwise alignments between the genome of interest (the reference) and all other leaves in the tree.

Also see the OSF Project page for the BRAG paper: https://osf.io/ak54t/

## Requirements

In addition to Python3, BRAG requires the following Python modules:

- numpy  
- pandas  
- scipy  
- statsmodels  
- ete3  
- matplotlib  
- biopython  

BRAG should be run on a machine with at least 16 GB of memory.

## Running BRAG with Docker

A docker image for BRAG is hosted on https://hub.docker.com/r/channsoden/brag/. You can use this image to run the sample analysis with:

```docker run channsoden/brag```

After the sample analysis has run you can copy the results out with:

```docker cp [container_name]:/home/BRAG/sample_data/results .```

Where [container_name] is the binomial tag generated for your container, viewed by:

```docker container ls -a```

Or you can use this image interactively to run your own analysis with:

```docker run -it channsoden/brag bash```

Consult the Docker documentation to mount filesystem directories on the container to get your data in/out: https://docs.docker.com/storage/bind-mounts/


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

This script will download (220 MB) and extract (710 MB) the genomes used in the BRAG paper, then recreate all figures and statistics from the paper. Running the analysis requires at least 8 GB memory.

