#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some phylogeny-related functions
(includeing one for solving half-breakpoint state parsimony problem)
"""

import math
from collections import defaultdict
from itertools import chain
import logging

import networkx as nx

from ragout.parsers.phylogeny_parser import (parse_tree, PhyloException)
from ragout.phylogeny.inferer import TreeInferer
from newick.tree import (Tree, Leaf)
logger = logging.getLogger()

class Phylogeny:
    """
    Represents phylogenetic tree and scores it with
    given half-breakpoint states
    """
    def __init__(self, tree):
        self.tree = tree
        self.tree_string = str(tree)
        self._scale_branches()

    @classmethod
    def from_newick(phylo, newick_str):
        return phylo(parse_tree(newick_str))

    @classmethod
    def from_permutations(phylo, perm_container):
        ti = TreeInferer(perm_container)
        return phylo(ti.build())

    def _scale_branches(self):
        """
        Fits mu coefficient according to branch lengths
        to avoid underflows/overflows
        """
        def tree_length(node):
            if node.terminal:
                return []

            lengths = []
            for node, _bootstrap, length in node.get_edges():
                assert length is not None
                lengths.append(length)
                lengths.extend(tree_length(node))

            return lengths

        lengths = tree_length(self.tree)
        assert len(lengths)
        self.mu = float(1) / _median(lengths)
        logger.debug("Branch lengths: {0}, mu = {1}".format(lengths, self.mu))

    def estimate_tree(self, leaf_states, internal_scores={}):
        """
        Scores the tree with weighted parsimony procedure
        """
        all_states = set(leaf_states.values())

        #score of a tree branch
        def branch_score(parent, child, branch):
            if parent == child or child is None:
                return 0.0
            else:
                #prevent underflow
                length = max(branch, 0.0000001)
                #adding one to counter possibly small exp value
                return 1.0 + math.exp(-self.mu * length)

        #recursive
        def rec_helper(root):
            if root.terminal:
                leaf_score = (lambda s: 0.0 if s == leaf_states[root.identifier]
                                            else float("inf"))
                return {s : leaf_score(s) for s in all_states}

            nodes_scores = {}
            for node, _bootstrap, _length  in root.get_edges():
                if node in internal_scores:
                    nodes_scores[node] = internal_scores[node]
                else:
                    nodes_scores[node] = rec_helper(node)

            root_scores = defaultdict(float)
            for root_state in all_states:
                for node, _bootstrap, branch_length in root.edges:
                    min_score = float("inf")
                    for child_state in all_states:
                        score = (nodes_scores[node][child_state] +
                                branch_score(root_state, child_state,
                                             branch_length))
                        min_score = min(min_score, score)
                    root_scores[root_state] += min_score

            return root_scores
        if self.tree in internal_scores.keys():
            return min(internal_scores[self.tree])
        return min(rec_helper(self.tree).values())

    def terminals_dfs_order(self):
        """
        Returns terminal nodes' names in dfs order
        """
        def get_labels(root):
            if root.terminal:
                return [root.identifier]

            edges = sorted(root.get_edges(), key=lambda e: e[2], reverse=True)
            edges = sorted(edges, key=lambda e: e[0].terminal)
            return list(chain(*map(lambda e: get_labels(e[0]), edges)))

        return get_labels(self.tree)

    def nodes_by_distance(self, genome, onlyLeaves=True):
        """
        Returns leaves names sorted by the distance from
        the given genome.
        """
        graph = nx.Graph()
        start = [None]
        def rec_helper(root):
            if root.identifier == genome:
                start[0] = root
            if root.terminal:
                return
            for node, _bootstrap, branch_length in root.edges:
                graph.add_edge(root, node, weight=branch_length)
                rec_helper(node)

        rec_helper(self.tree)
        distances = nx.single_source_dijkstra_path_length(graph, start[0])
        if onlyLeaves:
            nodes = [g for g in distances.keys()
                  if g.terminal and g.identifier != genome]
        else:
            nodes = [g for g in distances.keys()
                  if g.identifier != genome]
        return list(map(str, sorted(nodes, key=distances.get)))


def _median(values):
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) / 2]

def is_leaf(Node):
    if type(Node) == Leaf:
        return True
    return False

def get_node(Tree, identifier):
    if Tree.identifier != identifier:
        if is_leaf(Tree):
            return None
        else:
            for u, _bootstrap, length in Tree.get_edges():
                found = get_node(u, identifier)
                if found:
                    return found
                else:
                    continue
            return None
    else:
        return Tree

def estimate_labeled_tree(phylogeny, leaf_states, internal_states):
    """
    Scores the labeled-subtree with weighted parsimony procedure
    """
    all_states = set(leaf_states.values())

    #score of a tree branch
    def branch_score(parent, child, branch):
        if parent == child or child is None:
            return 0.0
        else:
            #prevent underflow
            length = max(branch, 0.0000001)
            #adding one to counter possibly small exp value
            return 1.0 + math.exp(-phylogeny.mu * length)

    #recursive
    def rec_helper(root):
        if root.terminal:
            leaf_score = (lambda s: 0.0 if s == leaf_states[root.identifier]
                                        else float("inf"))
            return {s : leaf_score(s) for s in all_states}

        nodes_scores = {}
        for node, _bootstrap, _length  in root.get_edges():
                nodes_scores[node] = rec_helper(node)

        root_scores = defaultdict(float)
        for root_state in all_states:

            if root.identifier in internal_states and root_state != internal_states[root.identifier]:
                root_scores[root_state] = float("inf")
                continue

            for node, _bootstrap, branch_length in root.edges:
                min_score = float("inf")
                for child_state in all_states:
                    score = (nodes_scores[node][child_state] +
                            branch_score(root_state, child_state,
                                         branch_length))
                    min_score = min(min_score, score)
                root_scores[root_state] += min_score
        return root_scores

    return rec_helper(phylogeny.tree)
