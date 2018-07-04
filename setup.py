from distutils.core import setup, Extension
from distutils.command.install import install as DistutilsInstall
import os

try:
    from Bio.PDB import *
except ImportError:
    print "\033[1;31mWARNING: Biopython installation not found\033[0m\n"

os.system('which g++ > .testgcc')
if os.stat(".testgcc").st_size == 0:
    print "\033[1;31mWARNING: g++ installation not found\033[0m\n"
os.remove(".testgcc")

class MyInstall(DistutilsInstall):
    def run(self):
        #do_pre_install_stuff()
        self.makeTefSearch()
        DistutilsInstall.run(self)
        #do_post_install_stuff()
    def makeTefSearch(self):
        os.chdir("src")
        os.system("make")
        os.chdir("..")

setup(
    cmdclass={'install': MyInstall},
    name = "TEF2",
    version = "2.0",
    #packages = ["tef2"],
    #scripts = ["tef2.py", "tef2/tef_search"],
    scripts = ['TEF2', 'src/tef_search'],
    packages = ["src"],
    #ext_modules=[Extension('tef_search', ['tef2/tef_search.cpp'])],
    description = "TEF",
    author = "Dirk Stratmann",
    author_email = "dirk.stratmann__at__upmc.fr",
    url = "http://www.impmc.upmc.fr",
    keywords = ["protein", "fragments", "closed loops", "tightened end fragments"],
    classifiers = [
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: Unix",
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    long_description = """\
TEF 2.0 

This tool decomposes proteins 3D structures into sub-domain fragments, called "tightened end fragments" (TEF) or "closed loops".	
(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J
Contact: dirk.stratmann@sorbonne-universite.fr


REQUIREMENTS
------------

Dependencies are:
- Biopython (http://www.biopython.org)
- (facultative) PyMOL
- (facultative) NACCESS (http://www.bioinf.manchester.ac.uk/naccess/)
"""
)

