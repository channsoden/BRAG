#!/usr/bin/env python
# Standard modules
import sys

# Nonstandard modules
import pandas as pd

# My modules
from fasta_tools import get_absolute_positions
from sakakibara_tools import read_sakakibara

def segment_tables(reference, segment_file, seqs_file, reference_genome_file):
    genomes, os_tab = read_sakakibara(segment_file, seqs_file, 'osfinder')

    # The ordering of the headers is given in seqs file, which contains the paths to the fasta file,
    # but the node labels in the tree should match the short names used in the filenames.
    # filenames are XXXXX.shortname1.shortname2.seqs.
    query = [name for name in seqs_file.split('.')[1:3] if name != reference]
    if len(query) != 1:
        sys.exit('Unable to identify reference in filename: {}'.format(seqs_file))
    else:
        query = query[0]
    pair = [line.strip() for line in open(seqs_file, 'r')] 
    ref_first = reference_genome_file in pair[0]
    ref_second = reference_genome_file in pair[1]
    if ref_first == ref_second:
        sys.exit('Unable to determine reference order in segment file {} with seqs file {}'.format(segment_file, seqs_file))

    if ref_first:
        osheader = ['rchr', 'rstart', 'rend', 'rsign', 'qchr', 'qstart', 'qend', 'qsign']
        reference_genome = genomes[0]
        query_genome = genomes[1]
    else:
        osheader = ['qchr', 'qstart', 'qend', 'qsign', 'rchr', 'rstart', 'rend', 'rsign']
        reference_genome = genomes[1]
        query_genome = genomes[0]
    os_tab.columns = osheader
    
    rscaffolds = scaffold_table(reference_genome)
    qscaffolds = scaffold_table(query_genome)

    add_absolute(os_tab, 'r', rscaffolds)
    add_absolute(os_tab, 'q', qscaffolds)

    #    os_tab['rcolors'] = alphabet[(os_tab.rchr-1) % len(alphabet)]
    #    os_tab['qcolors'] = alphabet[(os_tab.qchr-1) % len(alphabet)]

    # Signed indexes
    os_tab['ridx'] = os_tab.index.values * ((os_tab.rsign == '+') * 2 - 1)
    os_tab['qidx'] = os_tab.index.values * ((os_tab.qsign == '+') * 2 - 1)

    os_tab = os_tab.sort_values(['rstart_abs'])
    
    return query, rscaffolds, qscaffolds, os_tab

def scaffold_table(reference):
    scaffolds = get_absolute_positions(reference)
    scaffolds = {k.replace('NcraOR74A_', ''): v for k, v in list(scaffolds.items())}
    scaffolds = {'name': sorted(list(scaffolds.keys()), key=lambda k: scaffolds[k]), 'abs_pos': sorted(scaffolds.values())}
    scaffolds = pd.DataFrame.from_dict(scaffolds)
    return scaffolds

def add_absolute(os_tab, which, scaffolds):
    os_tab[which+'chr_start'] = scaffolds.iloc[(os_tab[which+'chr']-1).values].abs_pos.values
    os_tab[which+'start_abs'] = os_tab[which+'start'] + os_tab[which+'chr_start']
    os_tab[which+'end_abs'] = os_tab[which+'end'] + os_tab[which+'chr_start']
