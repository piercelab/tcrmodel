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

from TCR_functions import *

def create_profit_infile(infile):
    f = open(infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE *:*\n")
    f.write("FIT\n")
    f.close()

def get_rms(infile, reference, mobile):
    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    cmd = subprocess.Popen([profit_program, '-f', infile, '-h', reference, mobile], stdout=subprocess.PIPE)
    for line in cmd.stdout:
        rmsd = None
        if "RMS:" in line:
            rmsd = line.split(' ')[-1]
            return rmsd.rstrip()

def create_profit_cdr3_alignment_infile(profit_infile):
    if os.path.exists("fitted.pdb"):
        os.remove("fitted.pdb")
    f = open(profit_infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE 106-139:106-139\n")
    f.write("FIT\n")
    f.write("WRITE fitted.pdb\n")

def create_profit_cdr1_alignment_infile(profit_infile):
    if os.path.exists("fitted.pdb"):
        os.remove("fitted.pdb")
    f = open(profit_infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE 23-43:23-43\n")
    f.write("FIT\n")
    f.write("WRITE fitted.pdb\n")

def create_profit_cdr2hv4_alignment_infile(profit_infile):
    if os.path.exists("fitted.pdb"):
        os.remove("fitted.pdb")
    f = open(profit_infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE 55-91:55-91\n")
    f.write("FIT\n")
    f.write("WRITE fitted.pdb\n")

def create_profit_fw_alignment_infile(profit_infile):
    if os.path.exists("fitted.pdb"):
        os.remove("fitted.pdb")
    f = open(profit_infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE -23:-23\n")
    f.write("ZONE 43-55:43-55\n")
    f.write("ZONE 91-106:91-106\n")
    f.write("ZONE 139-:139-\n")
    f.write("FIT\n")
    f.write("WRITE fitted.pdb\n")
    
def align_using_cdr_segment(refpdb, mobpdb):
    profit_infile ="profit_cdr3.in"
    create_profit_cdr3_alignment_infile(profit_infile)
    rms = get_rms(profit_infile, refpdb, mobpdb)
    return rms

def find_cdr_template(cdrseq):
    identity_cutoff = 100
    scoring_matrix = "BLOSUM62"
    template_file = "/home/rg/benchmark/all_cdr3.seq"
    best_score = -9999999
    best_tag = ''
    best_seq =''
    with open(template_file, 'r') as f:
        for line in f:
            #print line
            linelist = line.split()
            if ( cdrseq == linelist[1] ): continue
            if ( len(cdrseq) == len(linelist[1]) ):
                if (identity_cutoff <= 100):
                    identity_score = calculate_identity_score( str(cdrseq), str(linelist[1]) )
                    #print "HI", linelist[2], linelist[3], identity_score
                    if ( calculate_identity_score(str(cdrseq), str(linelist[1])) > identity_cutoff ): continue
                    score = score_alignment(cdrseq,linelist[1],scoring_matrix)
                    if (score > best_score):
                        best_score = score
                        best_tag = linelist[0]
                        best_seq = linelist[1]
    #print best_score, best_tag
    mob = "/home/rg/tcr_template_update/update/" + best_tag + "_aho.pdb"
    ref = "/home/rg/tcr_template_update/update/" + native_tag + "_aho.pdb"
    #print ref, mob
    rms = align_using_cdr_segment(ref,mob)
    print cdrseq, len(cdrseq), best_score, best_tag, best_seq, rms

def create_profit_infile(infile):
    f = open(infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE *:*\n")
    f.write("FIT\n")
    f.close()

script, native_tag, cdrseq = argv
seqa_fw = ''
seqb_fw = ''

find_cdr_template(cdrseq)
