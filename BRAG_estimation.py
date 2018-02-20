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

class ref:
    # This is an empty container in the global namespace.
    # It is filled with data by set_reference().
    def __init__(self): pass

def set_reference(refnode, genome_length):
    ref.tree = refnode.get_tree_root()
    branches = list(ref.tree.iter_descendants())
    ref.tree_len = sum([n.dist for n in branches])
    ref.refnode = refnode
    ref.N = genome_length
    ref.masks = {branch:tree_mask(branch) for branch in branches}

class timer(object):
    def __init__(self):
        self.start = time.time()
        self.last = self.start

    def report(self):
        now = time.time()
        step = now - self.last
        total = now - self.start
        self.last = now
        return '{:.3f} seconds elapsed ({:.3f} total)'.format(step, total)
    
class tree_mask(object):
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
        
class qbreak(object):
    # A qbreak a region on the reference between two colinear sequences
    # that are adjacent in the reference and not adjacent in the query.
    # A qbreak is a vertex in the interval graph, and it's query is it's color.
    def __init__(self, query, start, end, certainty):
        self.query = query
        self.start = start
        self.end = end
        self.certain = certainty
        self.cliques = set()
        self.tbreaks = set()
        self.placed = False

    def __repr__(self):
        return '{}({},{})'.format(self.query, self.start, self.end)


    def place(self, tbreak):
        self.tbreaks.add(tbreak)

    def reset(self):
        assert not self.placed, 'should not reset placed qbreaks'
        self.tbreaks = set()
        self.cliques = set()

class clique(object):
    # A clique is a set of qbreaks that are all adjacent in the interval graph.
    # A clique is a contiguous region of the reference that has the same pattern of qbreaks,
    # i.e. a region of the reference where some qbreaks overlap.
    def __init__(self, start, end, qbreaks):
        self.start = start
        self.end = end
        self.qbreaks = qbreaks
        self.tbreaks = self.find_tbreaks()

        [qb.cliques.add(self) for qb in self.qbreaks]

    def __repr__(self):
        return '{}({},{})'.format(len(self.qbreaks), self.start, self.end)

    def find_tbreaks(self):
        # Returns most parsimonious locations of state transitions on tree (tbreaks)
        # given list of leaf nodes (self.qbreaks) that don't share the state of refnode,
        # assuming no convergent evolution to the reference state is possible.
        for i, qb in enumerate(self.qbreaks):
            if type(qb) == type(None):
                print 'this qb is a None:', i
                print len(self.qbreaks)
                print len([qb for qb in self.qbreaks if type(qb) == type(None)])
        different = set([ref.tree&qb.query for qb in self.qbreaks])
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
            evidence = [link_qbreaks(self.qbreaks, different - ori_leaves)]
        
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
            evidence.append(link_qbreaks(self.qbreaks, clade))
            different = different - clade
    
        assert set([qb for qbs in evidence for qb in qbs]) == set(self.qbreaks), 'not all qbreaks are used as evidence'
            
        tbreaks = set([tbreak(self, qbs, node) for qbs, node in zip(evidence, broken)])
        return tbreaks
    
    def edist(self):
        min_mask_lens = [ref.masks[tb.branch].minL for tb in self.tbreaks]
        max_mask_lens = [ref.masks[tb.branch].maxL for tb in self.tbreaks]
        min_edist = ref.tree_len - sum(max_mask_lens)
        max_edist = ref.tree_len - sum(min_mask_lens)
        return min_edist, max_edist

def link_qbreaks(qbreaks, leaves):
    # Link the broken branches to the qbreaks that support them
    leaf_names = [l.name for l in leaves]
    return [qb for qb in qbreaks if qb.query in leaf_names]
    
class tbreak(object):
    # Tree-consistent Break
    # A tbreak is a clique of qbreaks whose queries make up a monophyletic clade.
    # A tbreak's branch is the branch leading to the MRCA of the tbreak.
    def __init__(self, clique, qbreaks, branch):
        self.qbreaks = set(qbreaks)
        self.branch = branch
        self.cliques = set([clique])

        [qb.tbreaks.add(self) for qb in self.qbreaks]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            sort_key = lambda qb: qb.start
            return sorted(self.qbreaks, key=sort_key) == sorted(other.qbreaks, key=sort_key)
        return NotImplemented
    
    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def discard(self):
        [qb.tbreaks.remove(self) for qb in self.qbreaks]
        [clq.tbreaks.remove(self) for clq in self.cliques]

    def has_any_qb(self, qbreaks):
        for qb in self.qbreaks:
            if qb in qbreaks:
                return True
        return False

    def set_likelihood(self, block_length, size):
        self.likelihood = block_length / size


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

    tbreaks = merge_tbreaks(max_cliques)
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

    return rate_table


def merge_alignments(adj_coords):
    qbreaks = set([qbreak(query, *coord)
                   for query, coords in adj_coords
                   for coord in coords])
    return qbreaks

def find_maximal_cliques(qbreaks):
    # Find all maximal cliques in the interval graph of qbreaks.
    edges = [(ref.N, None)]
    for qb in qbreaks:
        edges.append((qb.start, qb))
        edges.append((qb.end, qb))
    edges.sort()

    last_removal = False
    open_qbreaks = set()
    cliques = [clique(0, edges[0][0], open_qbreaks)]
    for i, edge in enumerate(edges[:-1]):
        position, qb = edge
        clique_length = edges[i+1][0] - position
        try:
            open_qbreaks.remove(qb)
            if clique_length == 0 or last_removal:
                continue
            p = clique(position, position+clique_length, open_qbreaks)
            cliques.append(p)
            last_removal = True
        except KeyError:
            open_qbreaks.add(qb)
            last_removal = False

    cliques.append(clique(edges[-2][0], edges[-1][0], []))

    return cliques

def merge_tbreaks(cliques):
    # Find equivalent tbreaks and keep only one.
    tbreaks = cliques[0].tbreaks
    for i, clique in enumerate(cliques[:-1]):
        nxt = cliques[i+1]
        branches = set([tb.branch for tb in clique.tbreaks])
        nxt_branches = set([tb.branch for tb in nxt.tbreaks])
        duplicates = branches & nxt_branches
        new_tbs = [tb for tb in clique.tbreaks if tb.branch in duplicates]
        old_tbs = [tb for tb in nxt.tbreaks if tb.branch in duplicates]
        [old_tb.discard() for old_tb in old_tbs]
        tbreaks.update(nxt.tbreaks)
        nxt.tbreaks.update(new_tbs)
    return tbreaks

def find_subgraphs(tbreaks):
    # Find all connected subgraphs in the QT graph.
    qbs = set([max(tb.qbreaks, key=lambda qb: len(qb.tbreaks)) for tb in tbreaks])
    tb_subgraphs = []
    qb_subgraphs = []
    while qbs:
        qb = qbs.pop()
        sub_qbs, sub_tbs = follow_connections(qb)
        qbs -= sub_qbs
        tb_subgraphs.append(sub_tbs)
        qb_subgraphs.append(sub_qbs)
    return zip(tb_subgraphs, qb_subgraphs)

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
        qbs.update(*[tb.qbreaks for tb in unvisited])
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
            if tb1.qbreaks < tb2.qbreaks:
                subsets.add(tb1)
    [tb.discard() for tb in subsets]
    return tbreaks - subsets

def minimize_qbreaks(qbreaks):
    """Reduce the set of qbreaks within tbreaks to only subset-minimal qbreaks.
    A qbreak is subset-minimal if it's dual (qbreak.tbreaks) contains no other dual as a proper subset."""
    non_minimal = set()
    for qb1 in qbreaks:
        for qb2 in qbreaks:
            if qb1.tbreaks < qb2.tbreaks:
                non_minimal.add(qb2)
    # [tb.qbreaks.remove(qb) for qb in non_minimal for tb in qb.tbreaks]
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
            remaining_qbreaks = qbreaks - tb.qbreaks
            new_covered = covered_qbreaks | tb.qbreaks
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
    edges = [(ref.N, None, None)]
    for tb in tbreaks:
        start, end = place_tb(tb)
        edges.append((start, tb, None))
        edges.append((end, tb, None))
    for tb in tbs_in_sols.keys():
        start, end = place_tb(tb)
        edges.append((start, tb, tbs_in_sols[tb]))
        edges.append((end, tb, tbs_in_sols[tb]))
    edges.sort()
    
    open_tbreaks = []
    open_solutions = [[[] for cov in sol] for sol in solutions]
    # partiton = (start, end, length, [counts], [likelihoods], [tree_lengths])
    count_blocks = [(0, edges[0][0], edges[0][0], [0], [1], [ref.tree_len])]
    for i, edge in enumerate(edges[:-1]):
        position, tb, coord = edge
        end_pos = edges[i+1][0]
        block_length = float(end_pos - position)
        if coord:
            i, j = coord
            try:
                open_solutions[i][j].remove(tb)
            except ValueError:
                open_solutions[i][j].append(tb)
        else:
            try:
                open_tbreaks.remove(tb)
            except ValueError:
                open_tbreaks.append(tb)

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
            [tb.set_likelihood(block_length, tb.size) for tb in open_tbreaks]
            [tb.set_likelihood(block_length, tb.size) for sol in non_empty_solutions for cov in sol for tb in cov]

            cover_lengths = [[len(cov) for cov in sol] for sol in non_empty_solutions]
            num_placement_combos = sum([2 ** sum(combo) for combo in itertools.product(*cover_lengths)])
            if num_placement_combos > 10 ** 3:
                logfh.write('\tattempting to iterate over {} possible solutions'.format(num_placement_combos) )
            
            solution_combos = [x for x in itertools.product(*non_empty_solutions)]
            for s_combo in solution_combos:
                s_combo_tbreaks = open_tbreaks + [tb for cover in s_combo for tb in cover]
                
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

def place_tb(tb):
    start = max([qb.start for qb in tb.qbreaks])
    end = min([qb.end for qb in tb.qbreaks])
    tb.size = end - start
    return start, end

def estimate_rate(start, end, length, counts, likes, tree_lengths, rates):
    nucleotide_times = [length * tl for tl in tree_lengths]
    landscape = composite_poisson_likelihood(counts, likes, nucleotide_times, rates)
    return np.concatenate( (np.array([start, end, length]), landscape) )

def estimate_rates(count_blocks, rates, threads=1):
    estimations = [(estimate_rate, count_part + (rates,)) for count_part in count_blocks]
    estimates = mapPool(threads, estimations)
    header = ['start', 'end', 'length', 'E'] + [str(rate) for rate in rates]
    return pd.DataFrame(estimates, columns=header)    
