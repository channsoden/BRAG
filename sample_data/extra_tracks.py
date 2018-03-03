#!/usr/env/bin python
import sys
import numpy as np
import pandas as pd

from BRAG_parsers import scaffold_table

from fasta_tools import fasta_to_dict
from gff_tools import gff_table

def windowize(genome, step=4000, window_size=20000):
    start = 0
    end = start + window_size
    genlen = len(genome)
    while end < genlen:
        yield (start, end)
        start += step
        end += step

def flatten(df):
    # Merge overlapping regions into a single presence/absence layer.
    # Assumes df is sorted.
    regions = {'start':[], 'end':[]}
    region_start = 0
    region_end = 0
    for index, row in df.iterrows():
        if row.abs_start > region_end:
            regions['start'].append(region_start)
            regions['end'].append(region_end)
            region_start = row.abs_start
        region_end = row.abs_end
    regions['start'] = regions['start'][1:]
    regions['end'] = regions['end'][1:]
    return pd.DataFrame(regions)

def GC(seq):
    gc = seq.count('G') + seq.count('C')
    gc = float(gc) / len(seq)
    return gc

def CRI(seq):
    # composite rip index = (TA/AT)- (CA+TG)/(AC+GT)
    TA = float(seq.count('TA'))+1
    AT = float(seq.count('AT'))+1
    CA = float(seq.count('CA'))+1
    TG = float(seq.count('TG'))+1
    AC = float(seq.count('AC'))+1
    GT = float(seq.count('GT'))+1
    return (TA/AT) - (CA+TG)/(AC+GT)


def gene_density(annotation, window):
    start, end = window
    window_size = float(end - start)
    overlapping = annotation.loc[(annotation.start < end) & (annotation.end > start)]
    start = np.full( len(overlapping), start, dtype=int)
    starts = np.maximum(start, overlapping.start)
    end = np.full( len(overlapping), end, dtype=int)
    ends = np.minimum(end, overlapping.end)
    overlap = np.sum(ends - starts)
    return overlap / window_size

genome_file = 'genomes/Neurospora-crassa_OR74A_v12_fixed.fasta'
fasta = open(genome_file, 'r').read()
scaffolds = fasta[1:].split('>')
scaffolds = [s.split('\n', 1) for s in scaffolds]
scaffolds = [s[1].replace('\n', '') for s in scaffolds]
genome = ''.join(scaffolds)

annotation = gff_table('genomes/annotations/Neurospora_crassa_OR74A_v12_fixed.gtf')
scaffolds = scaffold_table('genomes/Neurospora-crassa_OR74A_v12_fixed.fasta')
abs_positions = {row['name']: row.abs_pos for i, row in scaffolds.iterrows()}
scaf_adjustment = annotation.seqname.replace(abs_positions)
annotation['abs_start'] = annotation.start + scaf_adjustment
annotation['abs_end'] = annotation.end + scaf_adjustment
annotation.sort_values('abs_start', inplace=True)

protein_coding = annotation.loc[annotation.feature == 'CDS']
protein_coding = flatten(protein_coding)
expressed = annotation.loc[annotation.feature == 'exon']
expressed = flatten(expressed)

columns = ['start', 'end', 'GC', 'CRI', 'cds_density', 'exon_density']
tracks = {column:[] for column in columns}
for window in windowize(genome):
    seg = genome[slice(*window)]
    seg = seg.upper()
    tracks['start'].append(window[0])
    tracks['end'].append(window[1])
    tracks['GC'].append( GC(seg) )
    tracks['CRI'].append( CRI(seg) )
    tracks['cds_density'].append( gene_density(protein_coding, window) )
    tracks['exon_density'].append( gene_density(expressed, window) )

tracks = pd.DataFrame(tracks)
tracks.to_csv('Ncra_extra_tracks.tsv', sep='\t', columns = columns, index=False)

# CRI and GC are colinear, basically the same thing, so both shouldn't be used for
# break rate regression. Exon density and CDS density are also highly related.
# Exon density is also heavily affected by RIP, so I want to include exon density
# only at un-ripped sites.
# Positive scores for CRI indicate RIP.
tracks['unripped_exon_density'] = tracks.exon_density
tracks.loc[tracks.CRI > 0, 'unripped_exon_density'] = np.NaN

tracks.to_csv('Ncra_extra_tracks_unique.tsv', sep='\t', columns = ['start', 'end', 'CRI', 'unripped_exon_density'], index=False)
