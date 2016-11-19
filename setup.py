import os
from setuptools import *

class BinaryMissingException(Exception):
    pass
if not os.path.isfile("lib/ragout-maf2synteny") or not os.path.isfile("lib/ragout-overlap"):
    raise BinaryMissingException("Maf2synteny bin doesn't exists, please run make first")
setup(name='RagoutAPI',
      version='0.1',
      description='API from Ragout for ancestor sequences reconstruction',
      url='http://github.com/ptdtan/Ragout',
      author='Ptdtan',
      author_email='ptdtan@gmail.com',
      license='MIT',
      packages=find_packages("."),
      install_requires = ['networkx'],
      scripts=['lib/ragout-maf2synteny', 'lib/ragout-overlap'],
      zip_safe=False)
