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
sys.path.insert(0,'/home/gowthamanr/tcr_update/template_update_scripts/')
from TCR_functions import *
script, infile = argv


print infile
cap_pdbfile = infile + "_tcrpmhc.pdb"
if os.path.isfile(cap_pdbfile): 
    pseq = pdbchain_to_fasta(cap_pdbfile, "C")
    plen= len(pseq)
    pfname = "MHC1a.P"+str(plen)+".fasta"
    pf=open(pfname, "a+")
    pf.write(">"+infile+"\n")
    mhcseq = pdbchain_to_fasta(cap_pdbfile, "A")
    pf.write(mhcseq+"\n")
