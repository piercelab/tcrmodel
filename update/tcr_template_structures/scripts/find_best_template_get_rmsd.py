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

    
script, native_tag, alphaseq, betaseq = argv

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

def create_profit_framework_alignment_infile(profit_infile):
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
    
def align_using_framework_segment(refpdb, mobpdb):
    profit_infile ="profit_framework.in"
    create_profit_framework_alignment_infile(profit_infile)
    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    myoutput = open("profit_outfile.out", 'w')
    processa = Popen([profit_program, '-f', profit_infile, '-h', refpdb, mobpdb], stdout=myoutput, stderr=PIPE)
    stdout, stderr = processa.communicate()
    #print stdout
    #print stderr

def find_orientation_template(seqa, seqb):
    identity_cutoff = 100
    scoring_matrix = "BLOSUM62"
    orientation_template_file = "/home/rg/tcr_template_update/update/TCR_ORIENTATION.seq"
    best_score = -9999999
    best_alpha_tag = ''
    best_beta_tag = ''
    #print seqa, seqb
    with open(orientation_template_file, 'r') as f:
        for line in f:
            #print line
            linelist = line.split()
            if ( (seqa == linelist[2]) and (seqb == linelist[3]) ): continue
            if ( (len(seqa) == len(linelist[2])) and (len(seqb) == len(linelist[3])) ):
                if (identity_cutoff <= 100):
                    identity_score = ( calculate_identity_score( str(seqa+seqb), str(linelist[2]+linelist[3]) ) )
                    #print "HI", linelist[2], linelist[3], identity_score
                    if (calculate_identity_score( str(seqa+seqb), str(linelist[2]+linelist[3]) ) > identity_cutoff ): continue
                    score = score_alignment(seqa,linelist[2],scoring_matrix) + score_alignment(seqb,linelist[3],scoring_matrix)
                    if (score > best_score):
                        best_score = score
                        best_alpha_tag = linelist[0]
                        best_beta_tag = linelist[1]
    print best_score, best_alpha_tag, best_beta_tag
    refa = "/home/rg/tcr_template_update/update/" + best_alpha_tag + "_aho.pdb"
    moba = "/home/rg/tcr_template_update/update/" + native_tag[:6]+ "_A_aho.pdb"
    align_using_framework_segment(refa,moba)
    copyfile("fitted.pdb", "fita.pdb")
    refb = "/home/rg/tcr_template_update/update/" + best_beta_tag + "_aho.pdb"
    mobb = "/home/rg/tcr_template_update/update/" + native_tag[:5] + native_tag[6:7] + "_B_aho.pdb"
    align_using_framework_segment(refb,mobb)
    copyfile("fitted.pdb", "fitb.pdb")
    os.system("cat fita.pdb fitb.pdb > fitab.pdb") 

def create_profit_infile(infile):
    f = open(infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE *:*\n")
    f.write("FIT\n")
    f.close()

seqa_fw = ''
seqb_fw = ''

print "START" , native_tag
if alphaseq:
    aseq = alphaseq    
    regexres = assign_CDRs_using_REGEX(aseq, 'A')
    if regexres:
        seqa_vdomain = regexres.groups()[0]
        seqa_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seqa_cdr1 = regexres.groups()[2]
        seqa_cdr2hv4 = regexres.groups()[4]
        seqa_cdr3 = regexres.groups()[6]   
        #print "cap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
    else:
        print "FAIL: RegEx failed to identify TCR domain A"


if betaseq:
    bseq = betaseq
    regexres = assign_CDRs_using_REGEX(bseq, 'B')
    if regexres:
        seqb_vdomain = regexres.groups()[0]
        seqb_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seqb_cdr1 = regexres.groups()[2]
        seqb_cdr2hv4 = regexres.groups()[4]
        seqb_cdr3 = regexres.groups()[6]   
        #print "cap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
    else:
        print "FAIL: RegEx failed to identify TCR domain B"

'''
print seqb_vdomain
print seqb_fw
print seqb_cdr1
print seqb_cdr2hv4
print seqb_cdr3
print seqa_fw,seqb_fw
'''

find_orientation_template(seqa_fw,seqb_fw)
infile = "infile.in"
create_profit_infile(infile)
native_pdb = "/home/rg/tcr_template_update/update/" + native_tag
rms = get_rms(infile, "fitab.pdb", native_pdb)
print native_tag, "rmsd", rms
