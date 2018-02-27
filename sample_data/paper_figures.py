#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy import stats
import statsmodels.stats.multitest as smm

from plots import regression_plot, pretty_bar
from plotting_tools import direct_labels
from BRAG_parsers import scaffold_table
from gff_tools import gff_table

def telomere_correlations(dataset, chromosomes, datalabel):
    dataset['midpoint'] = (dataset.start + dataset.end) / 2.
    tests = []
    for chr_num, coordinates in enumerate(chromosomes):
        label = 'Chr_{}_'.format(chr_num+1)
        lstart, lend, rstart, rend = coordinates # coordinates of left and right arms of the chromosomes

        larm = dataset.loc[(dataset.start >= lstart) & (dataset.end <= lend)]
        rarm = dataset.loc[(dataset.start >= rstart) & (dataset.end <= rend)]

        # left arm works from left to right
        ltest = regression_plot(larm.midpoint - lstart, larm[datalabel], label=label+'left')
        ltest.regress(slope = 'negative')
        tests.append(ltest)
        # right arm works from right to left
        rtest = regression_plot(rend - rarm.midpoint, rarm[datalabel], label=label+'right')
        rtest.regress(slope = 'negative')
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


certain = pd.read_csv('certain_rate_windows.txt', sep='\t', header=0, index_col=0)
uncertain = pd.read_csv('uncertain_rate_windows.txt', sep='\t', header=0, index_col=0)
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

# Conserved regions
uncertain.sort_values('start', inplace=True)
conserved = uncertain.loc[uncertain.E == 0]

low_regions = []
region_start = 0
region_end = 0
for index, row in conserved.iterrows():
    if row.start > region_end:
        low_regions.append( (region_start, region_end, region_end - region_start) )
        region_start = row.start
    region_end = row.end
low_regions = low_regions[1:]

lowfh = open('conserved_regions.txt', 'w')
conserved_length = 0
low_geneIDs = set()
for start, end, length in low_regions:
    conserved_length += int(length)
    features = annotation.loc[(annotation.abs_start <= end) & (annotation.abs_end >= start)]
    genes = set(features.attributes)
    low_geneIDs |= set([gene.split()[1][1:-2] for gene in genes])
    lowfh.write( 'Region {}-{} ({} bp): {} genes\n'.format(start, end, length, len(genes)) )
    for gene in genes:
        lowfh.write( '{}\n'.format(gene) )
    lowfh.write('\n')
lowfh.close()
print('{} genes found in {} bp of unbroken regions'.format(len(low_geneIDs), conserved_length))

# Most rapidly breaking regions
certain = certain.sort_values('E', ascending=False)
# get the top regions of the same length as the conserved regions
i = 0
high_length = 0
while high_length < conserved_length:
    i += 1
    high_windows = certain.head(i)
    high_length = np.sum( high_windows.end - high_windows.start )
high_windows = high_windows.sort_values('start')

high_regions = []
region_start = 0
region_end = 0
for index, row in high_windows.iterrows():
    if row.start > region_end:
        high_regions.append( (region_start, region_end, region_end - region_start) )
        region_start = row.start
    region_end = row.end
high_regions = high_regions[1:]

highfh = open('rapid_breakers.txt', 'w')
high_geneIDs = set()
for start, end, length in high_regions:
    features = annotation.loc[(annotation.abs_start <= end) & (annotation.abs_end >= start)]
    genes = set(features.attributes)
    high_geneIDs |= set([gene.split()[1][1:-2] for gene in genes])
    highfh.write( 'Region {}-{} ({} bp): {} genes\n'.format(start, end, length, len(genes)) )
    for gene in genes:
        highfh.write( '{}\n'.format(gene) )
    highfh.write('\n')
highfh.close()
print('{} genes found in {} bp of most rapidly breaking regions'.format(len(high_geneIDs), high_length))

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
lowcounts, label_order = core_enrichment(low_geneIDs, core_genes)
highcounts, label_order = core_enrichment(high_geneIDs, core_genes)
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
ax.set_xlim(-20, 140)
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

print('{} / {} genes in conserved regions have published phylogenetic distribution'.format(sum(lowcounts.values()), len(low_geneIDs)))
print('{} / {} genes in fragile regions have published phylogenetic distribution'.format(sum(highcounts.values()), len(high_geneIDs)))
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


# is break rate correlated with distance to telomere?

tests = telomere_correlations(certain, chromosomes, 'E')
#tests.extend( telomere_correlations(uncertain, chromosomes, 'E') )

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
        test.draw(ax, logy = logy, fit_report_location = (0.05, 0.05))
    fig.savefig(outfile, dpi=350)
    plt.close(fig)

grid_plots(tests, 'brag_correlations')


# is CRI correlated with distance to telomere?

print('regression of CRI w/ telomere distance')

extra_tracks = pd.read_csv('Ncra_extra_tracks.tsv', sep='\t', header=0)

tests = telomere_correlations(extra_tracks, chromosomes, 'CRI')

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

grid_plots(tests, 'CRI_correlations', logy=False)
