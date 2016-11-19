#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module implements a breakpoint graph
which is widely used in Ragout
"""

from itertools import chain
import os
import logging
from copy import copy
from collections import namedtuple

import networkx as nx
from networkx.drawing import nx_pydot as dot
from ragout.shared.debug import DebugConfig
from ragout.phylogeny.phylogeny import *

logger = logging.getLogger()
debugger = DebugConfig.get_instance()

GenChrPair = namedtuple("GenChrPair", ["genome", "chr"])

COLORS = ["blue", "green", "yellow", "black", "pink"]
COLORSDICT = dict()

class BreakpointGraph(object):
    """
    Breakpoint graph implementation
    """
    def __init__(self, perm_container=None, ancestral=False, ancestor=None, name="breakpoint_graph"):
        self.bp_graph = nx.MultiGraph()
        self.target = None
        self.references = []
        self.debug_nodes = set()
        self.ancestral = ancestral
        self.ancestor = ancestor
        if perm_container is not None:
            self.build_from(perm_container)
        self.writedot(name)

    def build_from(self, perm_container):
        """
        Builds breakpoint graph from permutations
        """
        for perm in perm_container.ref_perms:
            if perm.genome_name not in self.references:
                self.references.append(perm.genome_name)
        self.target = perm_container.target_perms[0].genome_name
        self.contig_ends = []
        for perm in perm_container.target_perms:
            self.contig_ends.append((perm.blocks[0].signed_id(),
                                     -perm.blocks[-1].signed_id()))

        c=0
        for perm in chain(perm_container.ref_perms,
                          perm_container.target_perms):
            assert perm.blocks
            if perm.genome_name not in COLORSDICT:
                COLORSDICT[perm.genome_name] = COLORS[c]
                c+=1
            for prev_block, next_block in perm.iter_pairs():
                self.bp_graph.add_node(-prev_block.signed_id())
                self.bp_graph.add_node(next_block.signed_id())

                self.bp_graph.add_edge(-prev_block.signed_id(),
                                       next_block.signed_id(),
                                       genome_id=perm.genome_name,
                                       chr_name=perm.chr_name,
                                       start=prev_block.end,
                                       end=next_block.start,
                                       infinity=False,
                                       color = COLORSDICT[perm.genome_name])

            if perm.genome_name in self.references and not perm.draft:
                self.bp_graph.add_edge(-perm.blocks[-1].signed_id(),
                                       perm.blocks[0].signed_id(),
                                       genome_id=perm.genome_name,
                                       chr_name=perm.chr_name,
                                       infinity=True,
                                       color = COLORSDICT[perm.genome_name])
            if perm.genome_name in self.target and not perm.draft and self.ancestral:
                self.bp_graph.add_edge(-perm.blocks[-1].signed_id(),
                                       perm.blocks[0].signed_id(),
                                       genome_id=perm.genome_name,
                                       chr_name=perm.chr_name,
                                       infinity=True,
                                       color = COLORSDICT[perm.genome_name])
            if self.target in perm.genome_name and not 'unlocalized' in perm.genome_name and not perm.draft and self.ancestral:
                self.bp_graph.add_edge(-perm.blocks[-1].signed_id(),
                                       perm.blocks[0].signed_id(),
                                       genome_id=perm.genome_name,
                                       chr_name=perm.chr_name,
                                       infinity=True,
                                       color = COLORS[c])
        if self.ancestral:
            for perm in perm_container.ancestor_perms:
                assert perm.blocks
                for prev_block, next_block in perm.iter_pairs():
                    self.bp_graph.add_node(-prev_block.signed_id())
                    self.bp_graph.add_node(next_block.signed_id())

                    self.bp_graph.add_edge(-prev_block.signed_id(),
                                       next_block.signed_id(),
                                       genome_id=perm.genome_name,
                                       chr_name=perm.chr_name,
                                       start=prev_block.end,
                                       end=next_block.start,
                                       infinity=False,
                                       color = "red")

        logger.debug("Built breakpoint graph with {0} nodes"
                                        .format(len(self.bp_graph)))

    def writedot(self, name):
        l = len(self.bp_graph.nodes())
        fout = open(name+".dot", "w")
        fout.write("graph G{\n")
        def _convert(n):
            if n>0:
                return '"%dh"' %abs(n)
            else:
                return '"%dt"' %abs(n)
        for n in self.bp_graph.nodes():
            fout.write("%s;\n" %(_convert(n)))
        a = set()
        for u,v in self.bp_graph.edges():
            if "".join(map(str,[u,v])) not in a:
                a.add("".join(map(str,[u,v])))
                for k,value in self.bp_graph[u][v].items():
                    fout.write("%s -- %s [color=%s, key=%d];\n" %(_convert(u),_convert(v),value["color"], k))
        fout.write("}\n")

    def connected_components(self):
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)
        bp_graphs = []
        for subgr in subgraphs:
            bg = BreakpointGraph()
            bg.target = self.target
            bg.references = copy(self.references)
            bg.ancestral =self.ancestral
            bg.ancestor = self.ancestor
            bg.bp_graph = subgr
            bp_graphs.append(bg)
        return bp_graphs

    def genomes_chrs_support(self, node_1, node_2):
        if not self.bp_graph.has_edge(node_1, node_2):
            return []
        return list(map(lambda e: GenChrPair(e["genome_id"], e["chr_name"]),
                    self.bp_graph[node_1][node_2].values()))

    def genomes_support(self, node_1, node_2):
        return list(map(lambda gp: gp.genome,
                    self.genomes_chrs_support(node_1, node_2)))

    def to_weighted_graph(self, phylogeny):
        """
        Converts a breakpoint graph into a weighted adjacency graph
        using half-breakpoint state parsimony problem
        """
        assert len(self.bp_graph) >= 2
        g = nx.Graph()
        g.add_nodes_from(self.bp_graph.nodes())

        for node in self.bp_graph.nodes():
            adjacencies = {}
            for neighbor in self.bp_graph.neighbors(node):
                for edge in self.bp_graph[node][neighbor].values():
                    adjacencies[edge["genome_id"]] = neighbor

            for ref_id in self.references:
                if ref_id not in adjacencies:
                    adjacencies[ref_id] = None  #"void" state in paper

            break_weights = {}
            for neighbor in self.bp_graph.neighbors(node):
                adjacencies[self.target] = neighbor
                break_weights[neighbor] = phylogeny.estimate_tree(adjacencies)

            #normalization
            total_weights = sum(break_weights.values())
            for neighbor in self.bp_graph.neighbors(node):
                weight = (break_weights[neighbor] / total_weights
                          if total_weights != 0 else 0)
                _update_edge(g, node, neighbor, weight)

        return g

    def to_weighted_graph2(self, phylogeny):
        """
        Converts a breakpoint graph into a weighted ancestral adjacency graph
        using half-breakpoint state parsimony problem
        """
        assert len(self.bp_graph) >= 2
        g = nx.Graph()
        g.add_nodes_from(self.bp_graph.nodes())
        #print self.references, self.target
        for node in self.bp_graph.nodes():
            adjacencies = {}
            for neighbor in self.bp_graph.neighbors(node):
                for edge in self.bp_graph[node][neighbor].values():
                    adjacencies[edge["genome_id"]] = neighbor

            for ref_id in self.references: #maybe unessesary
                if ref_id not in adjacencies:
                    adjacencies[ref_id] = None  #"void" state in paper
            if self.target not in adjacencies:
                adjacencies[self.target] = None #"void" state for target genome

            break_weights = {}
            for neighbor in self.bp_graph.neighbors(node):
                ancestor_state = {self.ancestor: neighbor} #assign state for ancestor
                ancestor_tree = get_node(phylogeny.tree, self.ancestor)
                ancestor_phylo = Phylogeny(ancestor_tree)
                ancestor_scores = estimate_labeled_tree(ancestor_phylo, adjacencies, ancestor_state)
                break_weights[neighbor] = phylogeny.estimate_tree(adjacencies, internal_scores={ancestor_tree: ancestor_scores})

            #normalization
            total_weights = sum(break_weights.values())
            for neighbor in self.bp_graph.neighbors(node):
                weight = (break_weights[neighbor] / total_weights
                          if total_weights != 0 else 0)
                _update_edge(g, node, neighbor, weight)

        return g
    """
    def add_debug_node(self, node):
        self.debug_nodes.add(node)
    """

    def alternating_cycle(self, node_1, node_2):
        """
        Determines if there is a cycle of alternating colors
        that goes through the given red-supported (!) edge
        """
        def get_genome_ids((u, v)):
            return self.genomes_support(u, v)

        good_path = False
        for path in self._alternating_paths(node_1, node_2):
            #assert len(path) % 2 == 0
            if len(path) % 2 != 0:
				#print len(path)
				return None
            if len(path) == 2:
                continue

            edges = list(zip(path[:-1], path[1:]))
            even_colors = list(map(get_genome_ids, edges[1::2]))
            even_good = all(map(lambda e: set(e) == set([self.target]),
                        even_colors))
            if not even_good:
                continue

            odd_colors = list(map(get_genome_ids, edges[0::2]))
            common_genomes = set(odd_colors[0])
            for edge_colors in odd_colors:
                common_genomes = common_genomes.intersection(edge_colors)

            if common_genomes:
                #self._check_distances(path)
                good_path = True
                break

        return len(path) / 2 if good_path else None

    """
    def _check_distances(self, path):
        assert len(path) % 2 == 0
        path.append(path[0])
        edges = list(zip(path[:-1], path[1:]))
        even_dist = list(map(lambda (n1, n2): self.get_distance(n1, n2),
                             edges[1::2]))
        odd_dist = list(map(lambda (n1, n2): self.get_distance(n1, n2),
                            edges[0::2]))
        diff = abs(sum(even_dist) - sum(odd_dist))
        coeff = float(diff) / (sum(even_dist) + sum(odd_dist))
        logger.debug(coeff)
    """

    def is_infinity(self, node_1, node_2):
        if not self.bp_graph.has_edge(node_1, node_2):
            return False

        for edge_data in self.bp_graph[node_1][node_2].values():
            if edge_data["infinity"]:
                return True
        return False

    def get_distance(self, node_1, node_2, phylogeny, ancestral=False):
        """
        Tries to guess the distance between synteny blocks
        in a target genome or ancestor genome
        """
        if not ancestral:
            genome = self.target
        else:
            genome = self.ancestor
        DEFAULT_DISTANCE = 0
        if not self.bp_graph.has_edge(node_1, node_2):
            return DEFAULT_DISTANCE
        distances = {e["genome_id"] : e["end"] - e["start"]
                     for e in self.bp_graph[node_1][node_2].values()}

        genomes_order = phylogeny.nodes_by_distance(genome, onlyLeaves = not ancestral )
        if len(distances.keys()) == 1 and genome == distances.keys()[0]:
            return distances[genome]
        for g in genomes_order:
            if g in distances:
                return distances[g]

    def debug_output(self):
        if not debugger.debugging:
            return

        graph_out = os.path.join(debugger.debug_dir, "breakpoint_graph.dot")
        _output_graph(self.bp_graph, graph_out)

    def _alternating_paths(self, src, dst):
        """
        Finds a path of alternating colors between two nodes
        """
        visited = set()
        def rec_helper(node, colored):
            if node == dst:
                return [[dst]]

            visited.add(node)
            paths = []
            for neighbor in self.bp_graph.neighbors(node):
                if neighbor in visited:
                    continue

                ##
                genomes = self.genomes_support(node, neighbor)
                non_target = set(filter(lambda g: g != self.target, genomes))
                if colored and len(non_target) == 0:
                    continue
                if not colored and self.target not in genomes:
                    continue
                ##

                far_paths = rec_helper(neighbor, not colored)
                map(lambda p: p.append(node), far_paths)
                paths.extend(far_paths)
            visited.remove(node)
            return paths

        paths = list(map(lambda p: p[::-1], rec_helper(src, True)))
        return paths


def _update_edge(graph, v1, v2, weight):
    """
    Helper function to update edge's weight
    """
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight


def _output_graph(graph, out_file):
    """
    Outputs graph in dot format
    """
    with open(out_file, "w") as fout:
        fout.write("graph {\n")
        for v1, v2, data in graph.edges_iter(data=True):
            fout.write("{0} -- {1}".format(v1, v2))
            if len(data):
                extra = list(map(lambda (k, v) : "{0}=\"{1}\"".format(k, v),
                                 data.items()))
                fout.write(" [" + ", ".join(extra) + "]")
            fout.write(";\n")
        fout.write("}")

