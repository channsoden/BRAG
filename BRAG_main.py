#!/usr/bin/env python
# Standard modules
import os, sys, argparse

# Nonstandard modules
import pandas as pd
import numpy as np
from scipy.stats import describe
from ete3 import Tree

# My modules
from processing_tools import mapPool
from tree_tools import root, tree_order
from BRAG_parsers import segment_tables
from BRAG_figures import degrading_coverage, OS_length_hist, plot_break_rate, correlation_scatter, track_correlation
import BRAG_estimation as br

def main(tree, reference, outgroup, output, segment_files, seqs_files,
         step=7000, window_size=35000, nthreads=1, centromeres=None, tracks=None):
    outfile = output+'.BRAG.stats'
    print 'Writing log to {}'.format(outfile)
    log = open(outfile, 'w')

    tree = Tree(tree)
    root(tree, [outgroup])

    reference_genome_file = infer_reference(seqs_files)
    table_jobs = [(segment_tables, (reference, segment_file, seqs_file, reference_genome_file))
                  for segment_file, seqs_file in zip(segment_files, seqs_files)]
    tables = mapPool(nthreads, table_jobs)

    order = tree_order(reference, tree)
    tables.sort(key=lambda t: order.index(t[0])) # sort by queries into order
    queries, rscaffolds, qscaffolds, os_tabs = zip(*tables)
    rscaffolds = rscaffolds[0] # all have same reference
    N = rscaffolds.iloc[-1].abs_pos # position of the end == reference genome size

    coverages = [np.sum(os_tab.rend - os_tab.rstart) / float(N) for os_tab in os_tabs]
    coverage_stats = describe(coverages)
    log.write('{} genomes aligned to {}.\n'.format(coverage_stats.nobs, reference))
    log.write('Minimum coverage:\t{}\n'.format(coverage_stats.minmax[0]))
    log.write('Mean coverage:\t{}\n'.format(coverage_stats.mean))
    log.write('Maximum coverage:\t{}\n'.format(coverage_stats.minmax[1]))
    log.write('SD coverage:\t{}\n'.format(coverage_stats.variance ** 0.5))
    log.write('Cumulative coverage:\t{}\n\n'.format(cumulative_coverage(os_tabs, N)))
    degrading_coverage(coverages, os_tabs, N, 'coverage_survival_curve.png')

    hist_jobs = [(OS_length_hist, (reference, query, os_tab)) for query, rscaffolds, qscaffolds, os_tab in tables]
    mapPool(nthreads, hist_jobs)

    log.write('\nEstimating break rates. . .\n\n')
    if not (os.path.isfile('uncertain_breakrate_rates.tab') and
            os.path.isfile('certain_breakrate_rates.tab') and
            os.path.isfile('uncertain_breakrate.log') and
            os.path.isfile('certain_breakrate.log')):
        adj_jobs = [(map_breakpoints, [os_tab]) for os_tab in os_tabs]
        uncertain_adj_coords = mapPool(nthreads, adj_jobs)
        certain_adj_coords = [[coord for coord in coords if coord[2]] for coords in uncertain_adj_coords]
        uncertain_adj_coords = zip(queries, uncertain_adj_coords)
        certain_adj_coords = zip(queries, certain_adj_coords)

        br.set_reference(tree&reference, N)

    log.write('Calculating Ambiguous (Unconfirmed) Break Rates:\n')
    if not (os.path.isfile('uncertain_breakrate_rates.tab') and
            os.path.isfile('uncertain_breakrate.log')):
        uncertain_estimates = br.break_rate(uncertain_adj_coords, output='uncertain_breakrate', threads=nthreads)
    else:
        uncertain_estimates = pd.read_csv('uncertain_breakrate_rates.tab', sep='\t')
    log.write(open('uncertain_breakrate.log', 'r').read())
    
    log.write('\nCalculating Unambiguous (Confirmed) Break Rates:\n')
    if not (os.path.isfile('certain_breakrate_rates.tab') and
            os.path.isfile('certain_breakrate.log')):
        certain_estimates = br.break_rate(certain_adj_coords, output='certain_breakrate', threads=nthreads)
    else:
        certain_estimates = pd.read_csv('certain_breakrate_rates.tab', sep='\t')
    log.write(open('certain_breakrate.log', 'r').read())

    if not os.path.isfile('uncertain_rate_windows.panda'):
        uncertain_rate_windows = rate_windows(uncertain_estimates, N, step=step, window_size=window_size)
        uncertain_rate_windows.to_csv('uncertain_rate_windows.txt', sep='\t', index=False)
    else:
        uncertain_rate_windows = pd.read_csv('uncertain_rate_windows.txt', sep='\t', header=0)

    if not os.path.isfile('certain_rate_windows.txt'):
        certain_rate_windows = rate_windows(certain_estimates, N, step=step, window_size=window_size)
        certain_rate_windows.to_csv('certain_rate_windows.txt', sep='\t', index=False)
    else:
        certain_rate_windows = pd.read_csv('certain_rate_windows.txt', sep='\t', header=0)

    # Mask Centromeres
    if centromeres:
        centromeres = [map( int, line.split('#')[0].strip().split() )
                       for line in open(centromeres, 'r')
                       if line.split('#')[0].strip()]
        abs_centromeres = [(rscaffolds.iloc[scaf_idx].abs_pos + start,
                            rscaffolds.iloc[scaf_idx].abs_pos + stop)
                           for scaf_idx, start, stop in centromeres]
        certain_rate_windows = mask(certain_rate_windows, abs_centromeres)
        uncertain_rate_windows = mask(uncertain_rate_windows, abs_centromeres)
    else:
        abs_centromeres = []

    # Mask scaffold edges
    scaffold_boundaries = [(x, x) for x in rscaffolds.abs_pos]
    certain_rate_windows = mask(certain_rate_windows, scaffold_boundaries, inclusive=False)
    uncertain_rate_windows = mask(uncertain_rate_windows, scaffold_boundaries, inclusive=False)

    # Add in extra data tracks and mask
    if tracks:
        tracks = pd.read_csv(tracks, sep='\t')
        tracks.sort_values('start', inplace=True)
        tracks = mask(tracks, scaffold_boundaries, inclusive=False)
    else:
        tracks = []
    track_labels = [label for label in list(tracks) if label not in ['start', 'end']]

    # Plot Figures
    print
    print 'Plotting break rates calculated with "True" qbreaks ("certain", lower bound estimate)'
    print 'against break rates calculated with "True" and "False" break rates ("uncertain:, upper bound estimate).'
    print output+'.uncertainty.png'
    indexer = certain_rate_windows['E'] != uncertain_rate_windows['E']
    correlation_scatter(certain_rate_windows['E'].loc[indexer],
                        uncertain_rate_windows['E'].loc[indexer],
                        output+'.uncertainty.png')

    if track_labels:
        print
        print 'Performing linear regression between extra data tracks and break rate.'
        print output+'.tracks_x_breakrate.png'
        track_correlation(certain_rate_windows, tracks, track_labels, output+'.tracks_x_breakrate.png')

    print
    print 'Plotting break rates and extra tracks along the reference genome.'
    print output + '.brMap.png'
    plot_break_rate(N, queries, os_tabs,
                    certain_estimates, uncertain_estimates,
                    certain_rate_windows, uncertain_rate_windows,
                    tracks, track_labels,
                    rscaffolds, abs_centromeres,
                    step, output + '.brMap.png')

def infer_reference(seqs_files):
    # The reference is the genome that appears in all the seqs files.
    # If there are multiple such genomes, or none, than something is wrong.
    potential_references = set([line.strip() for line in open(seqs_files[0], 'r')])
    for seqf in seqs_files:
        genomes = set([line.strip() for line in open(seqf, 'r')])
        potential_references = potential_references & genomes
    if len(potential_references) != 1:
        sys.exit('Reference genome could not be inferred from .seqs files provided. ' + \
                 'Either multiple genomes common in all alignments or no genome common to all alignments.')
    else:
        return list(potential_references)[0]

def cumulative_coverage(os_tabs, N):
    edges = []
    for tab in os_tabs:
        edges.extend(tab.rstart)
        edges.extend(-tab.rend)
    edges.sort(key=abs)

    covered = []
    icovered = False
    for edge in edges:
        if edge < 0 and icovered:
            end = abs(edge)
            covered.append((start, end))
            icovered = False
        elif edge >= 0 and not icovered:
            start = edge
            icovered = True
        else:
            pass
        
    coverage = sum([end-start for start, end in covered])
    return coverage / float(N)
    
def map_breakpoints(os_tab):
    rsort = os_tab.sort_values(['rstart_abs'])
    ref_adj, rheads, rtails = get_adjacencies('rchr', rsort)
    qsort = os_tab.sort_values(['qstart_abs'])
    qer_adj, qheads, qtails = get_adjacencies('qchr', qsort)

    # (start, end, True/False)
    # True coords are definitely breakpoints
    # False coords are possibly breakpoints, but unknown due to fragmentation of query assembly
    ref_adj_coords = [(os_tab.loc[os_tab.ridx == a, 'rend_abs'].iloc[0],
                    os_tab.loc[os_tab.ridx == b, 'rstart_abs'].iloc[0],
                    not ((a in qtails and b in qheads) or
                         (-a in qheads and -b in qtails)))
                   for a, b in ref_adj]
    
    return ref_adj_coords

def get_adjacencies(scaf_col, os_tab):
    adj = []
    heads = []
    tails = []
    scafs = set(os_tab[scaf_col])
    if scaf_col == 'rchr':
        idx = 'ridx'
    else:
        idx = 'qidx'
    for scaf in scafs:
        scaf = os_tab[os_tab[scaf_col] == scaf]
        heads.append(scaf[idx].iloc[0])
        tails.append(scaf[idx].iloc[-1])
        adj.extend(zip(scaf[idx][:-1], scaf[idx][1:]))
    return adj, heads, tails

def rate_windows(regions, N, window_size = 35000, step=7000):
    # Regions is a DataFrame w/ these columns:
    # 'start', 'end', 'length', 'E', rate_0. . . rate_50
    column_labels = regions.columns[5:]
    
    windows = pd.DataFrame(np.arange(0, N, step), columns=['start'])
    windows['end'] = windows.start + window_size

    def window_average(w):
        overlapping = regions.loc[(regions.start <= w.end) & (regions.end >= w.start)]
        missing_front = w.start - overlapping.start
        missing_front[missing_front < 0] = 0
        missing_back = overlapping.end - w.end
        missing_back[missing_back < 0] = 0
        overlap = overlapping.length - missing_front - missing_back
        window_expect = (overlap * overlapping['E']).sum() / window_size
        #columns = [(overlap * overlapping[label]).sum() / window_size for label in column_labels]
        #return tuple(columns)
        return window_expect

    windows['E'] = windows.apply(window_average, axis=1)

    return windows

def mask(windows, mask_intervals, inclusive=True):
    condition = np.zeros(len(windows))
    if inclusive:
        for start, end in mask_intervals:
            condition = np.logical_or(condition, np.logical_and(windows.start <= end, windows.end >= start))
    else:
        for start, end in mask_intervals:
            condition = np.logical_or(condition, np.logical_and(windows.start <= end, windows.end > start))
    return windows.mask(condition)


def genome_speeds(certain_rate_windows):
    # looks like maybe I'm trying to correlate rates with the distance to a centromere here?
    # I should finish this up, could be interesting
    # should look at each chromosome separately, as well as convert to % of total length and compare together
    certain_rate_windows['mean'] = (certain_rate_windows.min_rate + certain_rate_windows.max_rate) / 2

    certain_rate_windows['rscaf'] = 0
    for i, pos in zip(rscaffolds.abs_pos.index, rscaffolds.abs_pos):
        certain_rate_windows.loc[certain_rate_windows.start > pos, 'rscaf'] = i
    chromosomes = certain_rate_windows[np.in1d(certain_rate_windows.rscaf, [15, 6, 5, 3, 8, 12, 17])]
    chromosomes['middle'] = (chromosomes.start + chromosomes.end) / 2

    chromosomes['dist_to_cent'] = 0
    for scaf, start, end in centromeres:
        cent_middle = (end+start) / 2.
        chrom = chromosomes.loc[chromosomes.rscaf == scaf]
        chrom.dist_to_cent = (chrom.middle - cent_middle)
