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

from TCR_functions import *
script, pdb, aseq, bseq, date, resln, cdr3ab = argv



acdrfromseq = get_cdr_from_seq_by_aho_num(aseq,"A")
aseq_cdr3_aho = str(acdrfromseq[4])  

bcdrfromseq = get_cdr_from_seq_by_aho_num(bseq,"B")
bseq_cdr3_aho = str(bcdrfromseq[4])  

print pdb[:4], len(aseq_cdr3_aho), aseq_cdr3_aho, len(bseq_cdr3_aho), bseq_cdr3_aho
#print pdb, aseq, bseq, date, resln, aseq_cdr3_aho+bseq_cdr3_aho
