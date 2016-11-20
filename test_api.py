import os
from ragout import ragout_api as api

prefix = os.path.dirname(__file__)

ins = api.RagoutInstance(maf=os.path.join(prefix, "examples_API", "alignment.maf"),
                        ancestor="A",
                        ancestor_fasta = os.path.join(prefix, "examples_API", "ancestor.fasta"),
                        phylo="((Human:0.057804,Marmoset:0.04965)A:0.088104,B:0.064051)C;",
                        outDir=os.path.join(prefix, "test_API_output"),
                        scale="small",
                        is_overwrite=True)
ins._construct_ancestor()
