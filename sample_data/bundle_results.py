#!/usr/bin/env python

# Moves results to directory for easy export.
# Get container name with "sudo docker container ls -a"
# Copy out results with "sudo docker cp [NAME]:/home/BRAG/sample_data/results ."

import os
from glob import glob

# List input and configuration files.
exclude = ['genomes',
           'genomes.tar.gz',
           'alignments',
           'conserved_genes.py',
           'export_results.py',
           'extra_tracks.py',
           'matplotlibrc',
           'Neurospora-crassa_v12_fixed_centromeres.txt',
           'Neurospora.nwk',
           'paper_figures.py',
           'sample_analysis.sh',
           'takao_core_genes_key.txt',
           'takao_core_genes.tsv']

if not os.path.isdir('results'):
    os.mkdir('results')

files = [f[2:] for f in glob('./*')]
for f in files:
    if f not in exclude:
        os.rename(f, './results/'+f)
