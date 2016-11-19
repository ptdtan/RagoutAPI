#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module parses Ragout configuration file
"""

from collections import namedtuple
import re
import os
import logging

from ragout.parsers.phylogeny_parser import get_leaves_names, PhyloException
import ragout.shared.config as config

logger = logging.getLogger()

class RecipeException(Exception):
    pass

def _make_dummy_recipe(references, target, ancestor, tree, blocks, maf, naming_ref):
    lines = []
    lines.append("".join([".references = ", ",".join(references[1:])]))
    lines.append("".join([".tree = ", tree]))
    lines.append("".join([".target = ", target]))
    lines.append("".join([".ancestor = ", ancestor]))
    lines.append("".join([".blocks = ",blocks]))
    lines.append("".join([".maf = ", maf]))
    lines.append("".join([".naming_ref = ", naming_ref]))
    return parse_ragout_recipe(dummy_lines=lines)

def parse_ragout_recipe(filename=None, dummy_lines=None):
    if filename!=None:
        if not os.path.exists(filename):
            raise RecipeException("Can't open recipe file")

    if filename!=None:
        prefix = os.path.dirname(filename)
    else:
        prefix = ""

    recipe_dict = {"genomes" : {}}
    known_params = ["tree", "target", "blocks", "maf", "hal", "fasta",
                    "draft", "references", "naming_ref", "ancestor"]
    deprecated = ["circular"]
    required_params = ["references", "target"]

    cast_bool = ["circular", "draft"]
    fix_path = ["fasta", "maf", "hal"]

    defaults = {"circular" : False,
                "draft" : False}

    param_matcher = re.compile("([^\s]+)\s*=\s*([^\s].*)$")
    if dummy_lines:
        lines = dummy_lines
    else:
        lines = open(filename, "r").read().strip().split("\n")
    for lineno, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        m = param_matcher.match(line)
        if not m or not "." in m.group(1):
            raise RecipeException("Error parsing recipe on line {1}"
                                  .format(filename, lineno + 1))

        (obj, param_name), value = m.group(1).split("."), m.group(2)
        if param_name in deprecated:
            logger.warning("Recipe parameter '{0}' is deprecated"
                                                    .format(param_name))
            continue
        if param_name not in known_params:
            raise RecipeException("Unknown recipe parameter '{0}' on line {1}"
                                  .format(param_name, lineno, filename))

        #checking values, casting
        if param_name in cast_bool:
            if value.lower() in ["true", "1"]:
                value = True
            elif value.lower() in ["false", "0"]:
                value = False
            else:
                raise RecipeException("Error parsing recipe on line "
                                      "{0}: wrong value '{1}' for bool param"
                                      .format(lineno, value))
        if param_name == "blocks":
            if value not in config.vals["blocks"]:
                raise RecipeException("Unknown synteny block size set: {0}"
                                      .format(value))
        if param_name == "references":
            value = list(map(lambda s: s.strip(), value.split(",")))
        if param_name in fix_path:
            value = os.path.expanduser(value)
            try:
                value = os.path.join(prefix, value)
            except Exception as err:
                pass
        ###

        if obj == "":
            recipe_dict[param_name] = value
        elif obj == "*":
            defaults[param_name] = value
        else:
            recipe_dict["genomes"].setdefault(obj, {})[param_name] = value

    for param in required_params:
        if param not in recipe_dict:
            raise RecipeException("Required parameter '{0}' not found in recipe"
                                  .format(param))

    genomes = recipe_dict["references"] + [recipe_dict["target"]] + [recipe_dict["ancestor"]]
    if "tree" in recipe_dict:
        try:
            leaves = get_leaves_names(recipe_dict["tree"])
            if set(leaves) != set(genomes)-set([recipe_dict["ancestor"]]):
                raise RecipeException("The tree does not agree with "
                                      "the specified set of genomes")
        except PhyloException as e:
            raise RecipeException(e)

    for g in recipe_dict["genomes"]:
        if g not in genomes:
            raise RecipeException("Recipe error: genome '{0}' is not in "
                                  "specified as reference or target".format(g))

    for g in genomes:
        recipe_dict["genomes"].setdefault(g, {})

    for g, g_params in recipe_dict["genomes"].items():
        for def_key, def_val in defaults.items():
            g_params.setdefault(def_key, def_val)
    return recipe_dict
