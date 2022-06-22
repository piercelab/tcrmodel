#!/usr/bin/env python

import os

from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
from Bio import SeqIO
from TCR_functions import *

parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()


fasta_file =  "tcr_alpha.fasta"
#fasta_file =  "tcr_beta.fasta"

TAG = "A"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:   
    get_tcr_vdomain(record, TAG)
    
#get_pdb("1ao7", "AB", None, "/TCRmodeller/PDB_RELEASE/pdb_structures")

'''
#fasta_file =  "all_beta.fasta"
fasta_file =  "tcr_beta.fasta"
tag = "B"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    #tmp1a = assign_CDRs_using_REGEX(record.seq, tag)
    #if tmp1a is None:
        #print record.id, "FAIL"
        #print record.seq
        
    #print record.id, record.seq
    run_anarci(record.seq)
    #with open('anarci.out', 'r') as f:
        #lines=f.readlines()
    if os.path.isfile("anarci.out"):
        num_lines = sum(1 for line in open('anarci.out'))
        if num_lines > 5:
            with open('anarci.out', 'r') as f:
                lines=f.readlines()
                if lines[5].split('|')[2] is tag:
                    print record.id, lines[5].split('|')[2]
                
'''
