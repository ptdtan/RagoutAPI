import os
from ragout import ragout_api as api

prefix = os.path.dirname(__file__)

ins = api.RagoutInstance(maf=os.path.join(prefix, "examples_API", "alignment.maf"),
                        references=["Human", "Marmoset", "B"],
                        ancestor="A",
                        ancestor_fasta = os.path.join(prefix, "A", "ancestor.fasta"),
                        phyloStr="((Human:0.057804,Marmoset:0.04965)A:0.088104,B:0.064051)C;",
                        threads=4,
                        outDir=os.path.join(prefix, "test_API_output"),
                        scale="small",
                        is_overwrite=True)
ins._construct_ancestor()
