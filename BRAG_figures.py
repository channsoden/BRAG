#!/usr/bin/env python
# Standard modules

# Nonstandard modules
import numpy as np
import statsmodels.api as sm
import statsmodels.stats.multitest as smm
from matplotlib import pyplot as plt
from matplotlib import collections as mc
from matplotlib import patches, gridspec, colors
from matplotlib.colorbar import Colorbar

# My modules
from plotting_tools import alphabet, minimal, heat_scale, figure_add_alpha, simple_axis, percent_labels
from plots import regression_plot

### Histograms of orthologous segment lengths
def OS_length_hist(reference, query, os_tab):
    segment_lengths = os_tab.rend - os_tab.rstart
    inter_lengths = np.empty((0), dtype=int)
    scafs = set(os_tab.rchr)
    for scaf in scafs:
        scaf = os_tab[os_tab.rchr == scaf]
        inters = (scaf.rstart - scaf.rend.shift(1)).as_matrix()[1:]
        inter_lengths = np.concatenate((inter_lengths, inters), axis=0)

    hists = plt.figure(figsize=(8, 4))
    segL = hists.add_subplot(121)
    intL = hists.add_subplot(122)
    my_hist(segL, segment_lengths)
    my_hist(intL, inter_lengths)
    segL.set_xlabel('orthologous\nsegment length')
    intL.set_xlabel('non-orthologous\nsegment length')
    hists.savefig(reference+'_'+query+'_seg-lengths')
    plt.close(hists)

def my_hist(ax, series):
    bins = np.logspace(1, np.log10(series.max()), num=50)
    values, bins, patches = ax.hist(series, bins=bins, normed=True, color='grey', edgecolor='none')
    ax.set_xlim(series.min(), series.max())
    ax.set_xscale('log')

### Plot of global reference coverage by number of queries
def degrading_coverage(coverages, os_tabs, N, output):
    coverages, os_tabs = list(zip(*sorted(zip(coverages, os_tabs), reverse=True)))
    queries = list(range(len(os_tabs)))

    endlist = []
    for query, os_tab in enumerate(os_tabs):
        starts = [(start, query) for start in os_tab['rstart_abs']]
        ends = [(end, query) for end in os_tab['rend_abs']]
        endlist.extend(starts)
        endlist.extend(ends)
    endlist.sort()
    
    survivals = [N for q in queries]
    uncovered = set(queries)
    last_position = 0
    for position, query in endlist:
        if position != last_position:
            try:
                loser = min(uncovered)
                loss = position - last_position
                for q in queries[loser:]:
                    survivals[q] -= loss
            except ValueError:
                # all are covered, uncovered is empty
                pass
            last_position = position
        try:
            uncovered.remove(query)
        except KeyError:
            uncovered.add(query)

    survivals = [1]+[float(s)/N for s in survivals]

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(survivals, 'k-')

    ax.set_xlabel('Number of Queries')
    ax.set_ylabel('Coverage')
    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(queries))
    percent_labels(ax)
    simple_axis(ax)
    
    fig.savefig(output)
    plt.close(fig)
    
### Plot of rearrangement rates along reference genome
def plot_break_rate(N, queries, os_tabs,
                    certain_estimates, uncertain_estimates,
                    certain_rate_windows, uncertain_rate_windows,
                    extra_tracks, track_labels,
                    rscaffolds, abs_centromeres, step, outfile):
    nplots = 3 + len(track_labels)
    fig = plt.figure(figsize=(N/500000, 15))
    gs = gridspec.GridSpec(nplots, 1, height_ratios=[5, 3, 3] + [1]*len(track_labels))
    axes = [fig.add_subplot(gs[x]) for x in range(nplots)]
    ax = axes[0]
    certain_rate_ax = axes[1]
    uncertain_rate_ax = axes[2]

    [paint_reference(ax, i+1, os_tab) for i, os_tab in enumerate(os_tabs)]
    paint_breaks(ax, uncertain_rate_windows)
    rate_plot(certain_rate_ax, certain_estimates, certain_rate_windows, N, step=step, label='Certain\nBreak\nRate')
    rate_plot(uncertain_rate_ax, uncertain_estimates, uncertain_rate_windows, N, step=step, label='Uncertain\nBreak\nRate')

    window_x = (extra_tracks.start + extra_tracks.end) / 2
    for i, track in enumerate(track_labels):
        window_plot(axes[i+3], window_x, extra_tracks[track], track)
    
    # Chromosome Boundaries
    for i, scaf in rscaffolds.iterrows():
        [a.axvline(scaf.abs_pos, linewidth=2, linestyle='dashed') for a in axes]

    # Shade centromeres
    for abs_start, abs_stop in abs_centromeres:
        [a.axvspan(abs_start, abs_stop, color='black', alpha=0.5, lw=0) for a in axes]
    
    minimal(ax, labels=True)
    
    ax.set_yticks(list(range(0, len(os_tabs)+1)))
    ax.set_yticklabels([''] + list(queries))
    ax.set_ylim(-.6, len(os_tabs)+.6)
    
    xticks = (rscaffolds.abs_pos.iloc[:-1] + rscaffolds.abs_pos.shift(-1).iloc[:-1]) / 2
    for a in axes:
        a.set_xticks(xticks)
        a.set_xticklabels(rscaffolds.iloc[:-1].name)
        a.set_xlim(0, N)

    fig.savefig(outfile, dpi=350)
    plt.close(fig)
    
def paint_reference(ax, row, os_tab):
    #os_tab['rcolors'] = alphabet[(os_tab.rchr-1) % len(alphabet)]
    #os_tab['qcolors'] = alphabet[(os_tab.qchr-1) % len(alphabet)]
    
    segments = len(os_tab)
    lines = np.empty((segments, 2, 2), dtype=int)
    lines[:, :, 1] = row
    lines[:, 0, 0] = os_tab.rstart_abs
    lines[:, 1, 0] = os_tab.rend_abs
    
    lc = mc.LineCollection(lines, colors=alphabet[row % len(alphabet)], linewidths=10)
    ax.add_collection(lc)
    
def paint_breaks(ax, uncertain_rate_windows):
    # Heat map of the expected break rate.
    segments = len(uncertain_rate_windows)
    lines = np.empty((segments, 2, 2), dtype=int)
    lines[:, :, 1] = 0
    lines[:, 0, 0] = uncertain_rate_windows.start
    lines[:, 1, 0] = uncertain_rate_windows.end

    color_scale = uncertain_rate_windows['E'] / uncertain_rate_windows['E'].max()
    colors = heat_scale(color_scale)
    
    lc = mc.LineCollection(lines, colors=colors, linewidths=10)
    ax.add_collection(lc)

def rate_plot(ax, rate_estimates, rate_windows, N, step=7000, label=''):
    window_x = (rate_windows.start + rate_windows.end) / 2

    rates = [float(rate) for rate in rate_estimates.columns[5:]]
    log_rates = np.log10(rates)
    top = log_rates[-1]
    log_delta = log_rates[2] - log_rates[1]
    bottom = log_rates[1] - log_delta
    # enable drawing of estimates that overlap 0 (log10(0) = -inf)
    rates[0] = bottom
    log_expectation = np.log10(rate_windows.loc[:, 'E'])
    log_expectation[log_expectation < bottom] = bottom
    
    # Plot likelihood of rates as heat map
    extent = [0, N, bottom, top]
    likelihoods = rate_estimates.loc[:, rate_estimates.columns[5]:].T#.loc[::-1, :]
    hm_x = rate_estimates['start'].tolist() + [rate_estimates['end'].max()]
    heatmap = ax.pcolor(hm_x, log_rates, likelihoods, cmap='viridis',
                        norm=colors.LogNorm(vmin=1e-10, vmax=1))

    legend = label.replace('\n', '')+'_legend'
    legendfig = plt.figure(figsize=(.5,6))
    legend_ax = legendfig.add_subplot(111)
    cbar = Colorbar(legend_ax, heatmap)
    legendfig.savefig(legend)
    plt.close(legendfig)
    
    # Plot expectation as solid line
    ax.plot(window_x, log_expectation, 'k-', linewidth=2, label='Expected Rate')

    # legend = ax.legend(loc='upper left')
    ax.set_ylabel(label)
    ax.set_ylim(bottom, top)
    ax.yaxis.grid(which="major", color='k', linestyle='--', linewidth=1)

def window_plot(ax, x, y, label):
    ax.plot(x, y, 'k-', linewidth=2)
    ax.yaxis.grid(which="major", color='k', linestyle='--', linewidth=1)
    ax.set_ylabel(label)
    ax.set_ylim(np.min(y), np.max(y))

def correlation_scatter(x1, x2, outfile):
    y = x1 / x2

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    test = regression_plot(x2, y)
    test.regress()
    test.draw(ax, logx=True, fit_report_location = (0.05, 0.05))

    ax.set_ylabel('br(True) / br(True U False)')
    ax.set_xlabel('log( br(True U False) )')
    ax.set_ylim(0, 2)

    fig.savefig(outfile)
    plt.close(fig)
    return test

def track_correlation(rate_windows, tracks, track_labels, outfile):
    rate_x = (rate_windows.start + rate_windows.end) / 2
    track_x = (tracks.start + tracks.end) / 2
    assert np.sum(rate_x - track_x) == 0, 'Error: extra track window & step sizes do not match BRAG window & step sizes.'
    # Depending on how the extra tracks were windowized, there could be a different amount of
    # missing data windows (the meaningless end windows) than the BRAG rate_windows.
    # Since there is indexes match for non-missing data, reindex the shorter df with the longer's.
    if len(rate_windows) > len(tracks):
        tracks = tracks.reindex(rate_windows.index.values)
    if len(tracks) > len(rate_windows):
        rate_windodows = rate_windows.reindex(tracks.index.values)

    # Multiple linear regression.
    # I would think I might want to use a GLM, but I would have to know
    # something about the domain of the tracks. Since track values
    # could be any value from -inf to +inf, I think a Gaussian link function
    # might be the only thing suitable. Which I think is the same as MLR.
    # Additionally, the model coefficients will be unreliable when the tracks
    # are multicollinear.
    X = sm.add_constant(tracks[track_labels]) # allow for intercept
    model = sm.OLS(rate_windows.E, X, missing='drop')
    results = model.fit()
    print(results.summary())

    # Modeling each track with linear regression (simple, but most wrong).
    tests = []
    p_vals = []
    for label in track_labels:
        data = tracks[label]
        test = regression_plot(data, rate_windows.E, label = label)
        tests.append(test)
        p_vals.append( test.regress() )

    rejects, p_vals, bs, nonsense = smm.multipletests(p_vals, alpha=0.05, method= 'fdr_tsbh')
    for test, pv in zip(tests, p_vals):
        test.p_val = pv

    nplots = len(tests)
    columns = 4
    rows = int((nplots-1) / 2) + 1
    plotsize = (4, 4)
    figsize = (plotsize[0]*columns, plotsize[1]*rows)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(rows, columns)
    axes = [fig.add_subplot(gs[x]) for x in range(nplots)]

    for ax, test in zip(axes, tests):
        test.draw(ax, logy = True, fit_report_location = (0.05, 0.05), marker_alpha=0.1)

    fig.savefig(outfile)
    plt.close(fig)

