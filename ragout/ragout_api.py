
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

ragout_root = os.path.dirname(os.path.realpath(__file__))
lib_absolute = os.path.join(ragout_root, "lib")
sys.path.insert(0, lib_absolute)
sys.path.insert(0, ragout_root)
os.environ["PATH"] = lib_absolute + os.pathsep + os.environ["PATH"]

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
    Raogut instance for handling reconstruction methods
    """
    def __init__(self,maf, references, ancestor, ancestor_fasta,
                 threads=4,
                 phyloStr=None,
                 outDir="ragout-out",
                 scale="large",
                 tmpDir="tmp",
                 outLog="ragout-log.txt",
                 backend="maf",
                 is_overwrite = False,
                 is_debug=False,
                 is_resolve_repeats=False,
                 is_solid_scaffolds=False):
        self.maf = maf
        self.ancestor = ancestor
        self.ancestor_seqs = read_fasta_dict(ancestor_fasta)
        self.references = references
        self.target = references[0]
        self.phyloStr = phyloStr
        self.scale = scale
        self.debug = is_debug
        self.outDir = outDir
        self.backend = SyntenyBackend.backends[backend]
        self.overwrite = is_overwrite
        self.threads = threads
        if not tmpDir:
            self.tmpDir = os.path.join(outDir, "tmp")
        else:
            self.tmpDir = tmpDir

        self.phyloStr = phyloStr
        self.logger = enable_logging(outLog, is_debug)
        self.debugger = DebugConfig.get_instance()
        self.is_solid_scaffolds = is_solid_scaffolds
        self.is_resolve_repeats = is_resolve_repeats

        if not os.path.isdir(self.outDir):
            os.mkdir(self.outDir)
        if not os.path.isdir(self.tmpDir):
            os.mkdir(self.tmpDir)
        self.debug_root = self._set_debugging()
        self._set_exe_paths()
        self._check_extern_modules(backend)
        self.phylogeny, self.naming_ref = self._get_phylogeny_and_naming_ref()
        self.synteny_blocks = config.vals["blocks"][self.scale]
        self.dummy_recipe = _make_dummy_recipe(self.references, self.target, self.ancestor, self.phyloStr, self.scale, self.maf, self.naming_ref)
        self.perm_files = self._make_permutaion_files()
        self.run_stages = self.make_run_stages()
        self.phylo_perm_file = self.perm_files[self.synteny_blocks[-1]]
        self._make_stage_perms()
    def _construct_ancestor(self):

        ###Enable ChimeraDetector4Ancestor
        if not self.is_solid_scaffolds:
            raw_bp_graphs = {}
            for stage in self.run_stages:
                raw_bp_graphs[stage] = BreakpointGraph(self.stage_perms[stage],                          ancestor=self.ancestor, ancestral=True)
            chim_detect = ChimeraDetector4Ancestor(raw_bp_graphs, self.run_stages,  self.ancestor_seqs)

        prev_stages = []
        scaffolds = None
        ###apply for all stages
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

        ###output generating of ancestor scaffolds
        logger.info("Done scaffolding for ''{0}''".format(self.ancestor))
        out_gen = OutputGenerator(self.ancestor_seqs, scaffolds)
        out_gen.make_output(self.outDir, self.ancestor, write_fasta=False)

    def _set_debugging(self):
        if not os.path.isdir(self.outDir):
            os.mkdir(self.outDir)

        if not os.path.isdir(self.tmpDir):
            os.mkdir(self.tmpDir)

        debug_root = os.path.join(self.outDir, "debug")
        self.debugger.set_debugging(self.debug)
        self.debugger.set_debug_dir(debug_root)
        self.debugger.clear_debug_dir()
        return debug_root

    def _set_exe_paths(self, LIB_DIR="lib"):
        ragout_root = os.path.dirname(os.path.realpath(__file__))
        lib_absolute = os.path.join(ragout_root, LIB_DIR)
        sys.path.insert(0, lib_absolute)
        sys.path.insert(0, ragout_root)
        os.environ["PATH"] = lib_absolute + os.pathsep + os.environ["PATH"]
        pass

    def _check_extern_modules(self, backend):
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

    def _get_phylogeny_and_naming_ref(self):
        """
        Retrieves phylogeny (infers if necessary) as well as
        naming reference genome
        """
        if self.phyloStr:
            logger.info("Phylogeny is taken from parameters")
            phylogeny = Phylogeny.from_newick(self.phyloStr)
        else:
            raise Exception("Phylogeny tree must be supplied!")
            logger.info(phylogeny.tree_string)
        leaves_sorted = phylogeny.nodes_by_distance(self.target, onlyLeaves=True)
        naming_ref = leaves_sorted[0]
        logger.info("'{0}' is chosen as a naming reference".format(naming_ref))

        return phylogeny, naming_ref

    def _make_permutaion_files(self):
        return self.backend.make_permutations(self.dummy_recipe, self.synteny_blocks, self.outDir,
                                               self.overwrite, self.threads)

    def _make_stage_perms(self):
        self.stage_perms = {}
        for stage in self.run_stages:
            self.debugger.set_debug_dir(os.path.join(self.debug_root, stage.name))
            self.stage_perms[stage]= PermutationContainer(self.perm_files[stage.block_size],
                                                      self.dummy_recipe, stage.repeats,
                                                      stage.ref_indels, self.phylogeny)
            pass

    def make_run_stages(self):
        """
        Setting parameters of run stages
        """
        stages = []
        for block in self.synteny_blocks:
            stages.append(RunStage(name=str(block), block_size=block,
                                   ref_indels=False, repeats=False,
                                   rearrange=True))
        stages.append(RunStage(name="refine", block_size=self.synteny_blocks[-1],
                               ref_indels=False, repeats=self.is_resolve_repeats,
                               rearrange=False))
        return stages

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
