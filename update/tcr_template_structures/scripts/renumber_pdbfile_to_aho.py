#!/usr/bin/env python                                                                                                                                      
import sys, os, re
import gzip
from sys import argv
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()


import sys
sys.path.insert(0,'/home/rg/scripts')

import rosetta

from TCR_functions import *
script, pdb_file = argv
outfile = pdb_file + "aho.pdb"
renumber_pdbfile_to_aho(pdb_file,None,outfile,True)
