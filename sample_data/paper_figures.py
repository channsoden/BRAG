#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy import stats
# Suppress FutureWarning upon importing pandas.core.datetools in statsmodels
np.warnings.filterwarnings("ignore", category=FutureWarning)
import statsmodels.stats.multitest as smm

from plots import regression_plot, pretty_bar
from plotting_tools import direct_labels
from BRAG_parsers import scaffold_table
from gff_tools import gff_table

def telomere_correlations(dataset, chromosomes, datalabel, slope='negative'):
    dataset['midpoint'] = (dataset.start + dataset.end) / 2.
    tests = []
    for chr_num, coordinates in enumerate(chromosomes):
        label = 'Chr_{}_'.format(chr_num+1)
        lstart, lend, rstart, rend = coordinates # coordinates of left and right arms of the chromosomes

        larm = dataset.loc[(dataset.start >= lstart) & (dataset.end <= lend)]
        rarm = dataset.loc[(dataset.start >= rstart) & (dataset.end <= rend)]

        # left arm works from left to right
        ltest = regression_plot(larm.midpoint - lstart, larm[datalabel], label=label+'left')
        ltest.regress(slope = slope)
        tests.append(ltest)
        # right arm works from right to left
        rtest = regression_plot(rend - rarm.midpoint, rarm[datalabel], label=label+'right')
        rtest.regress(slope = slope)
        tests.append(rtest)
    return tests

def read_centromeres(cen_file, chromosomes):
    fh = open(cen_file, 'r')
    header = fh.readline()
    idxs, starts, stops = list(zip(*[list(map(int, line.split('#')[0].strip().split())) for line in fh]))
    centromeres = [[0, max(stops), 0, 0] for x in range(max(idxs)+1)]
    for idx, start, stop in zip(idxs, starts, stops):
        centromeres[idx][0] = chromosomes.iloc[idx].abs_pos
        centromeres[idx][1] = chromosomes.iloc[idx].abs_pos + min(start, centromeres[idx][1])
        centromeres[idx][2] = chromosomes.iloc[idx].abs_pos + max(stop, centromeres[idx][2])
        centromeres[idx][3] = chromosomes.iloc[idx+1].abs_pos -1
    return centromeres


certain = pd.read_csv('Ncrassa_certain_rate_windows.txt', sep='\t', header=0)
uncertain = pd.read_csv('Ncrassa_uncertain_rate_windows.txt', sep='\t', header=0)
# get the scaffold lengths
scaffolds = scaffold_table('genomes/Neurospora-crassa_OR74A_v12_fixed.fasta')
# calculate the coordinates of left and right arms of the chromosomes
chromosomes = read_centromeres('Neurospora-crassa_v12_fixed_centromeres.txt', scaffolds)

# look at gene content in high and low regions
abs_positions = {row['name']: row.abs_pos for i, row in scaffolds.iterrows()}
annotation = gff_table('genomes/annotations/Neurospora_crassa_OR74A_v12_fixed.gtf')
annotation = annotation.loc[annotation.feature == 'CDS'] # I only care about the CDS, not the start codons & exons etc
scaf_adjustment = annotation.seqname.replace(abs_positions)
annotation['abs_start'] = annotation.start + scaf_adjustment
annotation['abs_end'] = annotation.end + scaf_adjustment

def top_regions(data, annotation, top_size, outfile, ascending=False):
    data = data.sort_values('E', ascending=ascending)
    i = 0
    total_length = 0
    while total_length < top_size:
        i += 1
        top_windows = data.head(i)
        top_regions, total_length = collapse_regions(top_windows)

    fh = open(outfile, 'w')
    geneIDs = set()
    for start, end, length in top_regions:
        features = annotation.loc[(annotation.abs_start <= end) & (annotation.abs_end >= start)]
        genes = set(features.attributes)
        geneIDs |= set([gene.split()[1][1:-2] for gene in genes])
        fh.write( 'Region {}-{} ({} bp): {} genes\n'.format(start, end, length, len(genes)) )
        for gene in genes:
            fh.write( '{}\n'.format(gene) )
        fh.write('\n')
    fh.close()

    return geneIDs, total_length

def collapse_regions(data):
    data = data.sort_values('start')
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

most_size  = 2000000
# Conserved regions
low_genes, conserved_length = top_regions(certain, annotation, most_size, 'conserved_regions.txt', ascending=True)
print('{} genes found in {} bp of conserved regions'.format(len(low_genes), conserved_length))
high_genes, fragile_length = top_regions(certain, annotation, most_size, 'fragile_regions.txt', ascending=False)
print('{} genes found in {} bp of fragile regions'.format(len(high_genes), fragile_length))

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

    return category_counts, label_order

# Compare core gene content of high and low genes
print()
core_genes = parse_core_genes('takao_core_genes.tsv', 'takao_core_genes_key.txt')
lowcounts, label_order = core_enrichment(low_genes, core_genes)
highcounts, label_order = core_enrichment(high_genes, core_genes)
label_order = label_order[::-1]
diffcounts = [lowcounts[label]-highcounts[label] for label in label_order]


fig = plt.figure(figsize = (6, 6))
ax = fig.add_subplot(111)
bars = pretty_bar(ax, diffcounts, label_order, horizontal=True)
direct_labels(ax, list(range(len(diffcounts))), diffcounts, horizontal=True, num_format='0:+')
direct_labels(ax, list(range(len(diffcounts))),
              [((n<0)*2)-1 for n in diffcounts],
              altlabels=[highcounts[l] for l in label_order],
              horizontal=True)
ax.set_xlim(-105, 105)
ax.set_xlabel('Conserved Genes - Fragile Genes')
fig.savefig('core_gene_content')

core_genes_in_category = {}
for i, row in core_genes.iterrows():
    try:
        core_genes_in_category[row.ID30] += 1
    except KeyError:
        core_genes_in_category[row.ID30] = 1
print('{} genes have a published phylogenetic distribution'.format(sum(core_genes_in_category.values())))
for category, number in list(core_genes_in_category.items()):
    print('{}\t{}'.format(number, category))

print('{} / {} genes in conserved regions have published phylogenetic distribution'.format(sum(lowcounts.values()), len(low_genes)))
print('{} / {} genes in fragile regions have published phylogenetic distribution'.format(sum(highcounts.values()), len(high_genes)))
print()

distributions = [[lowcounts[label] for label in label_order],
                 [highcounts[label] for label in label_order]]
chi2, pval, dof, expected = stats.chi2_contingency(distributions)
print('2-sample Chi-squared test if the distribution of phylogenetic conservation of')
print('genes in fragile (most rapidly breaking) and conserved regions are equal:')
print('chi2', chi2)
print('pval', pval)
print('dof', dof)
print()

def multitest(df, chromosomes, track, slope='negative'):
    print('correlating {} with distance to telomere'.format(track))
    # is break rate correlated with distance to telomere?
    tests = telomere_correlations(df, chromosomes, track, slope=slope)
    
    # multiple testing correction for testing each chromosome arm independently
    # use 'fdr_tsbh' for two stage fdr correction via Benjamini/Hochberg method
    p_vals = [test.raw_slope_p for test in tests]
    rejects, p_vals, bs, nonsense = smm.multipletests(p_vals, alpha=0.05, method='fdr_tsbh')

    sig_tests = []
    insig_tests = []
    print('label\traw_p\tp\tslope\tintercept\tr_squared')
    for test, pv, reject in zip(tests, p_vals, rejects):
        test.p_val = pv
        print('{}\t{:.2E}\t{:.2E}\t{:.2E}\t{:.2E}\t{:0.4f}'.format(test.label, test.raw_slope_p, pv, test.slope, test.intercept, test.r2))
        if reject:
            sig_tests.append(test)
        else:
            insig_tests.append(test)

    print('sig tests', len(sig_tests))
    print('insig tests', len(insig_tests))

    return tests

def grid_plots(tests, outfile, logy=True):
    nplots = len(tests)
    columns = 2
    rows = int((nplots -1) / 2) + 1
    plotsize = (4, 4)
    figsize = (plotsize[0]*columns, plotsize[1]*rows)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(rows, columns)
    axes = [fig.add_subplot(gs[x]) for x in range(nplots)]

    for ax, test in zip(axes, tests):
        test.draw(ax, logy = logy, sig_color='g', fit_report_location = (0.05, 0.05))
    fig.savefig(outfile, dpi=350)
    plt.close(fig)

tests = multitest(certain, chromosomes, 'E')
grid_plots(tests, 'brag_correlations')

extra_tracks = pd.read_csv('Ncra_extra_tracks.tsv', sep='\t', header=0)
tests = multitest(extra_tracks, chromosomes, 'CRI')
grid_plots(tests, 'CRI_correlations', logy=False)
tests = multitest(extra_tracks, chromosomes, 'conservation', slope='positive')
grid_plots(tests, 'species_specific_genes_telomeres', logy=False)
