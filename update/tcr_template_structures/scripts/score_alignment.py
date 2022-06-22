#!/usr/bin/env python

import sys, os, re
import gzip
from sys import argv
import subprocess
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()
from shutil import copyfile
import glob

from TCR_functions import *

    
script, infile = argv


def get_best_align_score(inseq):
    identity_cutoff = 100
    scoring_matrix = "BLOSUM62"
    inlen = len(inseq)
    best_seq = ''
    best_score = -9999999
    best_tag = ''
    for name in glob.glob('/home/rg/tcr_template_update/update/*_CDR3_'+str(inlen)+'.fasta'):
        #print name
        for record in SeqIO.parse(name, "fasta"):
            if (identity_cutoff <= 100):
                identity_score = calculate_identity_score( str(record.seq), str(inseq) )
                if (identity_score > identity_cutoff ): continue
                score = score_alignment(str(record.seq), str(inseq), scoring_matrix)
                if (score > best_score):
                    best_seq = record.seq
                    best_score = score
                    best_tag = record.id
    print inseq, inlen, best_score, best_tag, best_seq

with open(infile, 'r') as f:
    for line in f:
        #strip newline
        newline = line.strip()
        #strip newline and first and last character
        #newline = line.strip()[1:-1]
        best_align_score = get_best_align_score(newline)
