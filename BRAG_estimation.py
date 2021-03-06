#!/usr/bin/env python
# Standard modules
import time, sys
import itertools

# Nonstandard modules
import numpy as np
import pandas as pd

# My modules
from BRAG_poisson import composite_poisson_likelihood
from processing_tools import mapPool

class timer:
    def __init__(self):
        self.start = time.time()
        self.last = self.start

    def report(self):
        now = time.time()
        step = now - self.last
        total = now - self.start
        self.last = now
        return '{:.3f} seconds elapsed ({:.3f} total)'.format(step, total)

class ref:
    # This is an empty container in the global namespace.
    # It is filled with data by set_reference().
    def __init__(self): pass

def set_reference(refnode, genome_length):
    ref.tree = refnode.get_tree_root()
    ref.branches = list(ref.tree.iter_descendants())
    ref.tree_len = 0
    i = 1
    for n in ref.branches:
        ref.tree_len += n.dist
        if not n.name:
            n.name = 'branch_{}'.format(i)
            i += 1
    ref.refnode = refnode
    ref.N = genome_length
    ref.masks = {branch:tree_mask(branch) for branch in ref.branches}
    
class tree_mask:
    def __init__(self, branch):
        descendants = list(branch.iter_descendants())
        if ref.refnode in descendants:
            # This branch is the origin of the reference state.
            masked = [node for node in ref.tree.traverse()
                      if node not in descendants]
        else:
            masked = descendants + [branch]
        max_mask_length = sum([n.dist for n in masked])
        # If the break is on the branch leading to the outgroup then the
        # root of the tree will bisect the branch to the outgroup.
        if branch.up.is_root():
            min_mask_length = max_mask_length - sum([n.dist for n in ref.tree.children])
        else:
            min_mask_length = max_mask_length - branch.dist

        self.leaves_masked = sum([n.is_leaf() for n in masked])
        self.minL = min_mask_length
        self.maxL = max_mask_length
        
class QBreak:
    # A qbreak a region on the reference between two colinear sequences
    # that are adjacent in the reference and not adjacent in the query.
    # A qbreak is a vertex in the interval graph, and it's query is it's color.

    # qbreak coordinates are [start, end) on nucleotide sequence
    # So broken bond could any between start and end sites
    def __init__(self, query, start, end, certainty):
        self.query = query
        self.start = start
        self.end = end
        self.certain = certainty
        self.cliques = set()
        self.tbreaks = set()

    def __hash__(self):
        return hash((self.query, self.start, self.end))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.query, self.start, self.end) == (other.query, other.start, other.end)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def __repr__(self):
        return '{}({},{})'.format(self.query, self.start, self.end)

class Clique(frozenset):
    # A clique is a set of qbreaks that are all adjacent in the interval graph.
    # A clique is a contiguous region of the reference that has the same pattern of qbreaks,
    # i.e. a region of the reference where some qbreaks overlap.
    def __new__(cls, qbreaks):
        return super().__new__(cls, qbreaks)

    def __init__(self, qbreaks):
        self.place()
        self.tbreaks = self.partition_tbreaks()
        [qb.cliques.add(self) for qb in self]

    def __repr__(self):
        return '{}({},{})'.format(len(self), self.start, self.end)

    def place(self):
        # Place the clique in the genome
        # Coordinates of qbreaks are [start, end)
        # So coordinates of a clique are also [start, end)
        self.start = max([qb.start for qb in self])
        self.end = min([qb.end for qb in self])
        self.size = self.end - self.start
        assert self.size >= 1, 'clique of size < 1 indicates invalid clique'

    def partition_tbreaks(self):
        # Finds most parsimonious locations of state transitions on tree (tbreaks)
        # given list of leaf nodes (self) that don't share the state of refnode,
        # assuming no convergent evolution to the reference state is possible.
        different = set([ref.tree&qb.query for qb in self])
        same = set(ref.tree.get_leaves()) - different # Leaves that share the state of refnode
        
        # Find the origin of the state of the reference
        origin = ref.refnode
        ori_leaves = origin.get_leaves()
        # while not all taxa sharing reference state are under the origin
        while not sum([taxon in ori_leaves for taxon in same]) == len(same):
            origin = origin.up
            ori_leaves = origin.get_leaves()
        ori_leaves = set(ori_leaves)
        
        if origin is ref.tree:
            # reference type is ancestral
            broken = []
            evidence = []
        else:
            broken = [origin]
            evidence = [link_qbreaks(self, different - ori_leaves)]
            
        # Find state transitions underneath the origin
        different = different & ori_leaves
        while different:
            node = next(iter(different))
            leaves = set(node.get_leaves())
            while leaves <= different:
                clade = leaves
                broken_node = node
                node = node.up
                try:
                    leaves = set(node.get_leaves())
                except AttributeError:
                    # node = None, we hit the root and tried to go up
                    break
            broken.append(broken_node)
            evidence.append(link_qbreaks(self, clade))
            different = different - clade
            
        assert set([qb for qbs in evidence for qb in qbs]) == self, 'not all qbreaks are used as evidence'
        
        tbreaks = [TBreak(qbs, node, self) for qbs, node in zip(evidence, broken)]
        return tbreaks

def link_qbreaks(qbreaks, leaves):
    # Link the broken branches to the qbreaks that support them
    leaf_names = [l.name for l in leaves]
    return [qb for qb in qbreaks if qb.query in leaf_names]

class TBreak(Clique):
    # Tree-consistent Break
    # A tbreak is a clique of qbreaks whose queries make up a monophyletic clade.
    # A tbreak's branch is the branch leading to the MRCA of the tbreak.
    def __new__(cls, qbreaks, branch, parent):
        return super().__new__(cls, qbreaks)

    def __init__(self, qbreaks, branch, parent):
        self.branch = branch
        self.clique = parent
        self.place()
        [qb.tbreaks.add(self) for qb in self]

    def __repr__(self):
        return '{}({},{})'.format(self.branch.name, self.start, self.end)

    def set_likelihood(self, block_length):
        self.likelihood = block_length / self.size

    def discard(self):
        # Remove the tbreak from it's qbreaks and parent clique
        # Needed to maintain duality when non-maximal tbreaks are tossed
        [qb.tbreaks.remove(self) for qb in self]
        self.clique.tbreaks.remove(self)

    def partition_tbreaks(self):
        # Override the clique method to avoid sensless recurrsion
        return None

def break_rate(adj_coords, output=None, threads=1):
    # Calculates the rearrangement rate at all break points
    for var in ['tree', 'tree_len', 'refnode', 'N', 'masks']:
        assert hasattr(ref, var), "Calculation of break rates requires '{}' to be set by set_reference()".format(var)

    if output:
        logfh = open(output+'.log', 'w')
    else:
        logfh = sys.stdout

    nodes = list(ref.tree.traverse())
    node_names = [str(i) for i in range(len(nodes))]
    leaves = [l for l in ref.tree.get_leaf_names()]

    logfh.write( 'Finding most parsimonious pattern of breaks. . .\n' )
    clock = timer()
    
    qbreaks = merge_alignments(adj_coords)
    logfh.write( '{}\tqbreaks\t{}\n'.format(len(qbreaks), clock.report()) )
    
    max_cliques = find_maximal_cliques(qbreaks)
    logfh.write( '{}\tmaximal_cliques\t{}\n'.format(len(max_cliques), clock.report()) )

    # Unique tbreaks are linked to qbreaks in sets.
    # tbreaks in cliques > tbreaks in qbreaks due to duplicates.
    tbreaks = frozenset([tb for qb in qbreaks for tb in qb.tbreaks])
    logfh.write( '{}\ttbreaks\t{}\n'.format(len(tbreaks), clock.report()) )
    
    sub_graphs = find_subgraphs(tbreaks)
    sub_graphs = [reduce_graph(*subgraph) for subgraph in sub_graphs]
    logfh.write( '{}\tdisconnected subgraphs\t{}\n'.format(len(sub_graphs), clock.report()) )
    
    essential_tbreaks = set()
    complex_sub_graphs = []
    for tbreaks, qbreaks in sub_graphs:
        reduced, simplicial, essential = remove_simplicial(qbreaks)
        essential_tbreaks.update(essential)
        if reduced:
            complex_sub_graphs.append((reduced, simplicial, essential))
    logfh.write( '{}\tsubgraphs with multiple solutions\t{}\n'.format(len(complex_sub_graphs), clock.report()) )

    solutions = [[cov - csg[2] for cov in minimum_covers(*csg)]
                 for csg in complex_sub_graphs]
    logfh.write( '{}\ttotal solutions found\t{}\n'.format(sum([len(sol) for sol in solutions]), clock.report()) )
    
    count_blocks = count_tbreaks(essential_tbreaks, solutions, logfh=logfh)
    logfh.write( '{}\tregions for which to estimate rates\t{}\n'.format(len(count_blocks), clock.report()) )

    header = ('start', 'end', 'length', 'counts', 'likelihoods', 'tree_lengths')
    count_table_file = output+'_counts.tab'
    count_table = pd.DataFrame(count_blocks, columns=header)
    count_table.to_csv(count_table_file, sep='\t')
    logfh.write( 'break counts written to {}\t{}\n'.format(count_table_file, clock.report()) )

    # maximum rate would be all leaves except reference and 1 other over 1 nucleotide
    max_breaks = len(ref.tree) - 2
    max_rate = max([ (max_breaks / ref.refnode.dist), (1. / (ref.refnode.dist / 2.)) ])
    # minimum rate would be 1 break over entire tree over entire genome
    min_rate = 1. / (ref.tree_len * ref.N)
    rates = np.concatenate( ([0], np.logspace(np.log10(min_rate), np.log10(max_rate))) )

    rate_table = estimate_rates(count_blocks, rates, threads=threads)
    logfh.write( 'break rates estimated for all regions\t{}\n'.format(clock.report()) )
    rate_table_file = output+'_rates.tab'
    rate_table.to_csv(rate_table_file, sep='\t')
    logfh.write( 'break rate estimates written to {}\t{}\n'.format(rate_table_file, clock.report()) )

    logfh.write( 'counting tbreaks by branch\t{}\n'.format(clock.report()) )
    tb_tree_file = output+'_tbreakTree'
    tb_tree = count_tbreaks_on_branches(essential_tbreaks, solutions, tb_tree_file)
    logfh.write( 'tree with tbreak counts as branch distances written to {}.nwk\n'.format(tb_tree_file) )
    logfh.write( 'tree with tbreak counts scaled by branch distances written to {}_scaled.nwk\t{}\n'.format(tb_tree_file, clock.report()) )

    return rate_table


def merge_alignments(adj_coords):
    qbreaks = set([QBreak(query, *coord)
                   for query, coords in adj_coords
                   for coord in coords])
    return qbreaks

def find_maximal_cliques(qbreaks):
    # Find all maximal cliques in the interval graph of qbreaks.
    # coordinates of qbreaks are [start, end)
    edges = []
    for qb in qbreaks:
        edges.append((qb.start, 1, qb))
        edges.append((qb.end, 0, qb))
    edges.sort(key=lambda edge: edge[0:2])

    open_qbreaks = set()
    cliques = []
    for i, edge in enumerate(edges):
        position, isstart, qb = edge
        if isstart:
            open_qbreaks.add(qb)
            next_edge = edges[i+1]
            if not next_edge[1]:
                # next edge is an end, so record this clique
                p = Clique(open_qbreaks)
                cliques.append(p)
        else:
            open_qbreaks.remove(qb)

    return cliques

def find_subgraphs(tbreaks):
    # Find all connected subgraphs in the QT graph.
    qbs = set([max(tb, key=lambda qb: len(qb.tbreaks)) for tb in tbreaks])
    tb_subgraphs = []
    qb_subgraphs = []
    while qbs:
        qb = qbs.pop()
        sub_qbs, sub_tbs = follow_connections(qb)
        qbs -= sub_qbs
        tb_subgraphs.append(sub_tbs)
        qb_subgraphs.append(sub_qbs)
    return list(zip(tb_subgraphs, qb_subgraphs))

def follow_connections(qbreak):
    qbs = set([qbreak])
    tbs = set()
    qbs_visited = set()
    tbs_visited = set()
    while qbs != qbs_visited:
        unvisited = qbs - qbs_visited
        tbs.update(*[qb.tbreaks for qb in unvisited])
        qbs_visited.update(unvisited)
        unvisited = tbs - tbs_visited
        qbs.update(*[tb for tb in unvisited])
        tbs_visited.update(unvisited)
    return qbs, tbs

def reduce_graph(tbreaks, qbreaks):
    # Do not change this order!
    return (maximize_tbreaks(tbreaks), minimize_qbreaks(qbreaks))

def maximize_tbreaks(tbreaks):
    """Remove tbreaks that are non-maximal.
    i.e. remove tbreaks that are subsets of other tbreaks."""
    subsets = set()
    for tb1 in tbreaks:
        for tb2 in tbreaks:
            if tb1 < tb2:
                subsets.add(tb1)
                tb1.discard()
                break
    return tbreaks - subsets

def minimize_qbreaks(qbreaks):
    """Reduce the set of qbreaks within tbreaks to only subset-minimal qbreaks.
    A qbreak is subset-minimal if it's dual (qbreak.tbreaks) contains no other dual as a proper subset."""
    non_minimal = set()
    for qb1 in qbreaks:
        for qb2 in qbreaks:
            if qb1.tbreaks < qb2.tbreaks:
                non_minimal.add(qb2)
    # [tb.remove(qb) for qb in non_minimal for tb in qb.tbreaks]
    return qbreaks - non_minimal
    
def remove_simplicial(qbreaks):
    """Find simplicial qbreaks, reduce the set of qbreaks to non-simplicial qbreaks,
    and identify the essential tbreaks necessitated by those qbreaks."""
    simplicial = set([qb for qb in qbreaks if len(qb.tbreaks) == 1])
    essential = set([tb for qb in simplicial for tb in qb.tbreaks])
    nonsimplicial = qbreaks - simplicial
    return nonsimplicial, simplicial, essential

def minimum_covers(qbreaks, covered_qbreaks, accepted_tbreaks):
    """Returns a list of all minimum tbreak covers of the qbreaks.
    Should be initialized with qbreaks being subset-minimal qbreaks,
    covered_qbreaks being the simplicial qbreaks,
    and accepted_tbreaks being the essential tbreaks.

    Complexity of this algorithm is about O(2**(len(tbreaks)-2)).
    May take an impossibly long time for certain unlikely worst-case scenarios."""
    if not qbreaks:
        return [accepted_tbreaks]
    else:
        covers = []
        qb = qbreaks.pop()
        for tb in qb.tbreaks:
            remaining_qbreaks = qbreaks - tb
            new_covered = covered_qbreaks | tb
            new_accepted = accepted_tbreaks | {tb}
            partial_covers = minimum_covers(remaining_qbreaks, new_covered, new_accepted)
            covers.extend(partial_covers)
        min_len = min([len(cov) for cov in covers])
        covers = [cover for cover in covers if len(cover) == min_len]
        return covers

def count_tbreaks(tbreaks, solutions, logfh=sys.stdout):
    # For each block formed by overlapping tbreaks:
    # Iterates through all possible tbreaks that could have occured in each region,
    # calculates the likelihood of each combination of tbreaks,
    # and returns the number of tbreaks and the observed tree length in each scenario,
    # as well as the likelihood of the scenario.

    # Each solution is a list of minimal covers for a disconnected subgraph.
    # While the set of minimal covers for each disconnected subgraph can be computed independently,
    # in the combined solution the disconnected subgraphs are not independent, since they could still
    # overlap spatially across the genome.
    tbs_in_sols = {tb:(i, j) for i, sol in enumerate(solutions) for j, cov in enumerate(sol) for tb in cov}

    # Find each block formed by overlapping tbreaks.
    # coordinates of tbreaks are [start, end)
    edges = [(ref.N, None, None)]
    for tb in tbreaks:
        edges.append((tb.start, tb, None))
        edges.append((tb.end, tb, None))
    for tb in list(tbs_in_sols.keys()):
        edges.append((tb.start, tb, tbs_in_sols[tb]))
        edges.append((tb.end, tb, tbs_in_sols[tb]))
    edges.sort(key = lambda edge: edge[0])
    
    open_tbreaks = set()
    open_solutions = [[set() for cov in sol] for sol in solutions]
    # partiton = (start, end, length, [counts], [likelihoods], [tree_lengths])
    count_blocks = [(0, edges[0][0], edges[0][0], [0], [1], [ref.tree_len])]
    for i, edge in enumerate(edges[:-1]):
        position, tb, coord = edge
        next_edge = edges[i+1]
        end_pos = next_edge[0]
        block_length = float(end_pos - position)
        if coord:
            i, j = coord
            try:
                open_solutions[i][j].remove(tb)
            except KeyError:
                open_solutions[i][j].add(tb)
        else:
            try:
                open_tbreaks.remove(tb)
            except KeyError:
                open_tbreaks.add(tb)

        if block_length == 0:
            continue

        # the break rate per base per unit time within this block
        # breaks per base per substitution
        min_rates = []
        max_rates = []

        counts = []
        likes = []
        tree_lengths = []
        
        non_empty_solutions = [sol for sol in open_solutions if sum([len(cov) for cov in sol])]
        
        if open_tbreaks or non_empty_solutions:
            # set likelihood of tbreaks occuring within this block
            [tb.set_likelihood(block_length) for tb in open_tbreaks]
            [tb.set_likelihood(block_length) for sol in non_empty_solutions for cov in sol for tb in cov]

            cover_lengths = [[len(cov) for cov in sol] for sol in non_empty_solutions]
            num_placement_combos = sum([2 ** sum(combo) for combo in itertools.product(*cover_lengths)])
            if num_placement_combos > 10 ** 3:
                logfh.write('\tattempting to iterate over {} possible solutions'.format(num_placement_combos) )
            
            solution_combos = [x for x in itertools.product(*non_empty_solutions)]
            for s_combo in solution_combos:
                open_combos = set([tb for cover in s_combo for tb in cover])
                s_combo_tbreaks = open_tbreaks | open_combos
                
                # likelihood = P(these solutions are true) * P(the breaks occured in this block, rather than another location)
                # likelihood = P(combo | solutions) * P(placements | tbreaks)
                certain = set([tb for tb in s_combo_tbreaks if tb.likelihood == 1])
                uncertain = set([tb for tb in s_combo_tbreaks if tb.likelihood != 1])

                # Find the power set of tbreaks that are not certain to be in this region
                placement_combos = [set(p_combo) for r in range(len(uncertain)+1)
                                    for p_combo in itertools.combinations(uncertain, r)]
                for p_combo in placement_combos:
                    num_breaks = len(certain) + len(p_combo)
                    
                    complement = uncertain - p_combo
                    placement_likelihood = np.product([tb.likelihood for tb in p_combo]) * np.product([1.-tb.likelihood for tb in complement])
                    composite_likelihood = placement_likelihood * (1. / len(solution_combos))
                    
                    min_observed_tree_length = max_observed_tree_length = ref.tree_len
                    for tb in p_combo | certain:
                        min_observed_tree_length -= ref.masks[tb.branch].maxL
                        max_observed_tree_length -= ref.masks[tb.branch].minL
                    
                    counts.append(num_breaks)
                    likes.append(composite_likelihood)
                    # Use the mean of the min and max observed tree lengths
                    # since given that one event occured in the interval of the branch length
                    # the probability density function of the event occuring in any part of the branch
                    # is uniform, assuming this is a poisson process
                    # Therefore, the expected time of occurance on the branch is the midpoint.
                    tree_lengths.append( (min_observed_tree_length + max_observed_tree_length) / 2. )
                
        else:
            # no breaks observed in this region
            counts.append(0)
            likes.append(1.)
            tree_lengths.append(ref.tree_len)
        
        part = (position, end_pos, block_length, counts, likes, tree_lengths)
        count_blocks.append(part)
        
    return count_blocks

def count_tbreaks_on_branches(tbreaks, solutions, outfile):
    tbob = {br.name:0 for br in ref.branches}
    for tb in tbreaks:
        tbob[tb.branch.name] += 1

    for sol in solutions:
        prob = 1. / len(sol)
        for cover in sol:
            for tb in cover:
                tbob[tb.branch.name] += prob

    tb_tree = ref.tree.copy()
    for branch in tb_tree.iter_descendants():
        branch.dist = tbob[branch.name]
    tb_tree_scaled = ref.tree.copy()
    for branch in tb_tree_scaled.iter_descendants():
        branch.dist = tbob[branch.name] / branch.dist

    tb_tree.write(outfile=outfile+'.nwk')
    tb_tree_scaled.write(outfile=outfile+'_scaled.nwk')
    return tb_tree

def estimate_rate(start, end, length, counts, likes, tree_lengths, rates):
    nucleotide_times = [length * tl for tl in tree_lengths]
    landscape = composite_poisson_likelihood(counts, likes, nucleotide_times, rates)
    return np.concatenate( (np.array([start, end, length]), landscape) )

def estimate_rates(count_blocks, rates, threads=1):
    estimations = [(estimate_rate, count_part + (rates,)) for count_part in count_blocks]
    estimates = mapPool(threads, estimations)
    header = ['start', 'end', 'length', 'E'] + [str(rate) for rate in rates]
    return pd.DataFrame(estimates, columns=header)    
