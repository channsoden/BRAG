#!/usr/bin/env python
import sys, os, random
sys.path.append('/home/christopher/brag_dev/BRAG/')
sys.path.append('/home/christopher/brag_dev/BRAG/hannsoden-bioinformatics/')
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import collections as mc
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy import stats
# Suppress FutureWarning upon importing pandas.core.datetools in statsmodels
np.warnings.filterwarnings("ignore", category=FutureWarning)
import statsmodels.stats.multitest as smm

from plots import regression_plot, pretty_bar
from plotting_tools import direct_labels, alphabet, brewer
from BRAG_parsers import scaffold_table
from gff_tools import gff_table

certain = pd.read_csv('Ncrassa_certain_rate_windows.txt', sep='\t', header=0)
uncertain = pd.read_csv('Ncrassa_uncertain_rate_windows.txt', sep='\t', header=0)
# get the scaffold lengths

# calculate the coordinates of left and right arms of the chromosomes
scaffolds = scaffold_table('genomes/Neurospora-crassa_OR74A_v12_fixed.fasta')
genome_size = int(scaffolds.loc[scaffolds.name == 'end', 'abs_pos'])
# look at gene content in high and low regions
abs_positions = {row['name']: row.abs_pos for i, row in scaffolds.iterrows()}
annotation = gff_table('genomes/annotations/Neurospora_crassa_OR74A_v12_fixed.gtf')
annotation = annotation.loc[annotation.feature == 'CDS'] # I only care about the CDS, not the start codons & exons etc
scaf_adjustment = annotation.seqname.replace(abs_positions)
annotation['abs_start'] = annotation.start + scaf_adjustment
annotation['abs_end'] = annotation.end + scaf_adjustment

def top_regions(data_by_E, data_by_start, annotation, top_size, startat=0, ascending=False):
    data = data_by_E
    i = startat
    total_length = 0
    while total_length < top_size:
        i += 1
        top_windows = data_by_start[data_by_start['start'].isin(data.head(i)['start'])]
        top_regions, total_length = collapse_regions(top_windows)

    geneIDs = set()
    for start, end, length in top_regions:
        features = annotation.loc[(annotation.abs_start <= end) & (annotation.abs_end >= start)]
        genes = set(features.attributes)
        geneIDs |= set([gene.split()[1][1:-2] for gene in genes])

    return geneIDs, total_length, i

def collapse_regions(data):
    regions = []
    total_length = 0
    region_start = 0
    region_end = 0
    for index, row in data.iterrows():
        if row.start > region_end:
            region_length = region_end - region_start
            total_length += region_length
            regions.append( (region_start, region_end, region_length) )
            region_start = row.start
        region_end = row.end
    regions = regions[1:]
    return regions, total_length

def parse_core_genes(gene_table, key_file):
    gene_table = pd.read_csv(gene_table, delimiter='\t', header=0)
    keyfh = open(key_file, 'r')
    ranks = {}
    better_labels = {}
    for line in keyfh:
        rank, code, description = line.strip().split('\t')
        ranks[code] = int(rank)
        better_labels[code] = description
    gene_table['coreness'] = gene_table.ID30
    gene_table.coreness.replace(ranks, inplace=True)
    gene_table.ID30.replace(better_labels, inplace=True)
    gene_table.sort_values('coreness', inplace=True)
    gene_table = gene_table.loc[gene_table.coreness != 6] # Drop "other" rank
    
    return gene_table

def core_enrichment(genes, core_genes):
    core_genes_present = core_genes.loc[core_genes.Broad7_geneID.isin(genes)]

    category_counts = {}
    label_order = []
    for i, row in core_genes_present.iterrows():
        try:
            category_counts[row.ID30] += 1
        except KeyError:
            category_counts[row.ID30] = 1
            label_order.append(row.ID30)
    label_order = [l for l in label_order if l != 'Others']

    return category_counts

core_genes = parse_core_genes('takao_core_genes.tsv', 'takao_core_genes_key.txt')

# Windows sorted from fragile to conserved
certain.sort_values('E', ascending=False, inplace=True)
categories = ['N. crassa', 'Ascomycota', 'Pezizomycotina', 'Dikarya', 'Cellular Life']
#colors = ['red', 'orange', 'yellow', 'green', 'blue']
#colors = [(1-(i/4.), 0, (i/4.)) for i in range(5)]
#colors = brewer[:5]
#colors = [alphabet[i] for i in [0, 1, 2, 3, 5]]
colors = ['#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6'][::-1]

if True:
    kw_ranks = {cat:[] for cat in categories}
    #rates = np.linspace(certain.E.iloc[-1], certain.E.iloc[0], num=50)
    rates = np.logspace(-4.5, -2, num=50, base=10)
    rates = np.append(rates[::-1], [-1])
    rate_hists = {cat:[] for cat in categories}

    total_genes = set()
    current_rank = 0
    rate_idx = 1
    current_rate = rates[rate_idx]
    new_hist_content = {cat:0 for cat in categories}
    for i, row in certain.iterrows():
        features = annotation.loc[(annotation.abs_start <= row.end) & (annotation.abs_end >= row.start)]
        attributes = set(features.attributes)
        genes = set([attr.split()[1][1:-2] for attr in attributes])
        new_genes = genes - total_genes
        total_genes |= new_genes
        new_content = core_enrichment(new_genes, core_genes)

        rate = row.E
        while rate <= current_rate:
            for cat, count in new_hist_content.items():
                rate_hists[cat].append( count )
            rate_idx += 1
            current_rate = rates[rate_idx]
            new_hist_content = {cat:0 for cat in categories}
        
        tie_value = (current_rank + len(new_genes)) / 2.
        current_rank += len(new_genes)
        for cat, count in new_content.items():
            kw_ranks[cat].extend( [tie_value] * count )
            new_hist_content[cat] += count
            
    for cat, count in new_hist_content.items():
        rate_hists[cat].append( count )

data = [kw_ranks[cat] for cat in categories]
fig = plt.figure(figsize = (6, 6))
ax = fig.add_subplot(111)
values, bins, patches = ax.hist(data, 30, histtype='bar', stacked=True, label=categories, color=colors)
label_x = 2400
last_mean = 0
for vals, label in zip(values, categories):
    mean = np.mean(vals)
    label_y = (mean + last_mean) / 2
    last_mean = mean
    ax.text(label_x, label_y, label, ha='center', va='bottom')
#ax.legend()
ax.set_ylabel('Genes')
ax.set_xlabel('Break Rate Rank')
ax.invert_xaxis()
fig.savefig('core_genes_2.png', dpi=350)

print('Kruskal-Wallis test for equal distributions')
kw_result = stats.kruskal(*kw_ranks.values())
print(kw_result)
for category in categories:
    print('{} mean rank:\t{}'.format(category, np.mean(kw_ranks[category])))

