
"""
Ragout API interface for ancestral genome reconstruction
"""
import os
import sys
import shutil
import logging
import argparse
from collections import namedtuple
from copy import deepcopy

import ragout.assembly_graph.assembly_refine as asref
import ragout.scaffolder.scaffolder as scfldr
import ragout.scaffolder.merge_iters as merge
import ragout.maf2synteny.maf2synteny as m2s
import ragout.overlap.overlap as overlap
import ragout.shared.config as config
from ragout.scaffolder.output_generator import OutputGenerator
from ragout.overlap.overlap import OverlapException
from ragout.phylogeny.phylogeny import Phylogeny, PhyloException
from ragout.breakpoint_graph.permutation import (PermutationContainer,
                                                 PermException)
from ragout.synteny_backend.synteny_backend import (SyntenyBackend,
                                                    BackendException)
from ragout.parsers.recipe_parser import parse_ragout_recipe, RecipeException, _make_dummy_recipe
from ragout.parsers.fasta_parser import read_fasta_dict, FastaError
from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import (Permutation, Block, Contig, Scaffold, Link)
from ragout.breakpoint_graph.breakpoint_graph import BreakpointGraph
from ragout.breakpoint_graph.inferer import AdjacencyInferer
from ragout.breakpoint_graph.chimera_detector import ChimeraDetector
from ragout.breakpoint_graph.chimera_detector_ancestor import ChimeraDetector4Ancestor
from ragout.phylogeny.phylogeny import *
from ragout.__version__ import __version__

import ragout.synteny_backend.maf


RunStage = namedtuple("RunStage", ["name", "block_size", "ref_indels",
                                   "repeats", "rearrange"])

class RagoutInstance(object):
    """
    Ragout instance for handling reconstruction methods
    """
    def __init__(self, maf, ancestor, ancestor_fasta, phylo,
                 outDir="ragout_out",
                 scale="large",
                 outLog="ragout_log",
                 is_overwrite = True,
                 is_debug = True,
                 is_resolve_repeats=False,
                 is_solid_scaffolds=False):
        self.maf = maf
        self.ancestor = ancestor
        self.ancestor_seqs = read_fasta_dict(ancestor_fasta)
        self.scale = scale
        self.debug = is_debug
        self.outDir = "_".join([outDir, ancestor])
        self.backend = SyntenyBackend.backends["maf"]
        self.overwrite = is_overwrite
        self.logger = enable_logging("_".join([outLog, ancestor]), is_debug)
        self.debugger = DebugConfig.get_instance()
        self.is_solid_scaffolds = is_solid_scaffolds
        self.is_resolve_repeats = is_resolve_repeats
        self.synteny_blocks = config.vals["blocks"][self.scale]

        if not os.path.isdir(self.outDir):
            os.mkdir(self.outDir)

        self._set_debugging()
        self._check_extern_modules()
        self._get_phylogeny_and_naming_ref(phylo)
        self._make_recipe(phylo)
        self._make_permutaion_files()
        self._make_run_stages()
        self._make_stage_perms()
        pass

    def _construct_ancestor(self):

        ###Enable ChimeraDetector4Ancestor
        if not self.is_solid_scaffolds:
            raw_bp_graphs = {}
            for stage in self.run_stages:
                raw_bp_graphs[stage] = BreakpointGraph(self.stage_perms[stage],                          ancestor=self.ancestor, ancestral=True)
            chim_detect = ChimeraDetector4Ancestor(raw_bp_graphs, self.run_stages,  self.ancestor_seqs)

        prev_stages = []
        scaffolds = None

        ###Iterative scaffolding
        last_stage = self.run_stages[-1]
        for stage in self.run_stages:
            logger.info("Stage \"{0}\"".format(stage.name))
            #debugger.set_debug_dir(os.path.join(debug_root, stage.name))
            prev_stages.append(stage)

            if not self.is_solid_scaffolds:
                broken_perms = chim_detect.break_contigs(self.stage_perms[stage], [stage])
            else:
                broken_perms = self.stage_perms[stage]
            breakpoint_graph = BreakpointGraph(broken_perms, ancestral=True, ancestor=self.ancestor)
            adj_inferer = AdjacencyInferer(breakpoint_graph, self.phylogeny, ancestral= True)
            adjacencies = adj_inferer.infer_adjacencies()
            cur_scaffolds = scfldr.build_scaffolds(adjacencies, broken_perms, ancestral=True)

            if scaffolds is not None:
                if not self.is_solid_scaffolds:
                    merging_perms = chim_detect.break_contigs(self.stage_perms[stage],
                                                              prev_stages)
                else:
                    merging_perms = self.stage_perms[stage]
                scaffolds = merge.merge_scaffolds(scaffolds, cur_scaffolds,
                                                  merging_perms, stage.rearrange, ancestral=True)
            else:
                scaffolds = cur_scaffolds
        scfldr.assign_scaffold_names(scaffolds, self.stage_perms[last_stage], self.naming_ref)

        ###Output generating of ancestor scaffolds
        logger.info("Done scaffolding for ''{0}''".format(self.ancestor))
        out_gen = OutputGenerator(self.ancestor_seqs, scaffolds)
        out_gen.make_output(self.outDir, self.ancestor, write_fasta=False)
        pass

    def _set_debugging(self):
        """
        Set debugging config
        """
        self.debug_root = os.path.join(self.outDir, "_   debug")
        self.debugger.set_debugging(self.debug)
        self.debugger.set_debug_dir(self.debug_root)
        self.debugger.clear_debug_dir()
        pass

    def _check_extern_modules(self):
        """
        Checks if all necessary native modules are available
        """
        if not m2s.check_binary():
            raise BackendException("maf2synteny binary is missing, "
                                   "did you run 'make'?")
        if not overlap.check_binary():
            raise BackendException("overlap binary is missing, "
                                   "did you run 'make'?")
        pass

    def _get_phylogeny_and_naming_ref(self, phyloStr):
        """
        Retrieves phylogeny (infers if necessary) as well as
        naming reference genome
        """
        if phyloStr:
            logger.info("Phylogeny is taken from parameters")
            self.phylogeny = Phylogeny.from_newick(phyloStr)
            if len(self.phylogeny.tree.get_leaves_identifiers()) !=3:
                raise Exception("The API only support for 3-leaves phylogeny tree")
        else:
            raise Exception("Phylogeny tree must be supplied!")
        logger.info(self.phylogeny.tree_string)
        self.references = self.phylogeny.tree.get_leaves_identifiers()
        self.target = self.references[0]
        leaves_sorted = self.phylogeny.nodes_by_distance(self.target, onlyLeaves=True)
        self.naming_ref = leaves_sorted[0]
        logger.info("'{0}' is chosen as a naming reference".format(self.naming_ref))
        pass

    def _make_permutaion_files(self):
        """
        Run maf2synteny to make permutation files from maf alignment
        """
        self.perm_files = self.backend.make_permutations(self.dummy_recipe, self.synteny_blocks, self.outDir, self.overwrite)
        pass

    def _make_stage_perms(self):
        self.stage_perms = {}
        for stage in self.run_stages:
            self.debugger.set_debug_dir(os.path.join(self.debug_root, stage.name))
            self.stage_perms[stage]= PermutationContainer(self.perm_files[stage.block_size],
                                                      self.dummy_recipe, stage.repeats,
                                                      stage.ref_indels, self.phylogeny)
            pass

    def _make_run_stages(self):
        """
        Setting parameters of run stages
        """
        self.run_stages = []
        for block in self.synteny_blocks:
            self.run_stages.append(RunStage(name=str(block), block_size=block,
                                   ref_indels=False, repeats=False,
                                   rearrange=True))
        self.run_stages.append(RunStage(name="refine", block_size=self.synteny_blocks[-1],
                               ref_indels=False, repeats=self.is_resolve_repeats,
                               rearrange=False))
        pass

    def _make_recipe(self, phylo):
        """
        Make a dummy recipe to fit Ragout logic
        """
        self.dummy_recipe = _make_dummy_recipe(self.references, self.target, self.ancestor, phylo, self.scale, self.maf, self.naming_ref)
        pass

def enable_logging(log_file, debug):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    logger = logging.getLogger()
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)
    return logger
