#!/usr/bin/env python

import sys, os, re
import gzip
from sys import argv
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
from Bio.PDB.Polypeptide import three_to_one

parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()

def run_anarci(seq, anarci_prog):
    myseq = str(seq.rstrip())
    process = Popen(["ANARCI", "-i", myseq, "-o", "anr_seq.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    anarcilist = []

    regex0 = "\#\|.*\|([A-B])\|.*\|.*\|([0-9]+)\|([0-9]+)\|"

    TAG = "~!@#$%^&*()_+"#just random initialization
    with open("anr_seq.out") as f:
        for line in f:
            res = re.search(regex0, line)
            if res:
                tcrdomain = myseq[int(res.groups()[1]):int(res.groups()[2])+1]
                if len(tcrdomain) == len(myseq):
                    TAG = res.groups()[0]

            if line.startswith(TAG):
                values = [line[0],line[2:7],line[8],line[10]]
                if line[10] == "-":continue
                anarcilist.append(values)

    return anarcilist

#####
#MAIN
#####
script, infile = argv
anarci_path = 'ANARCI'
inpstruct = parser.get_structure('TCR', infile)
ahostruct = inpstruct
    
#ppb = PPBuilder()# Using C-N
ppb = CaPPBuilder()# Using CA-CA
for chaininfo in inpstruct.get_chains():
    s = inpstruct[0][chaininfo.id]
    pp = ppb.build_peptides(s)
    seq = pp[0].get_sequence()
    anrlist = run_anarci(seq, anarci_path)
    i=0
    for residue in s:
        assert (three_to_one(residue.resname) == anrlist[i][3]), "Residues does not match %r %r %r" %(i, three_to_one(residue.resname), anrlist[i][3])
        residue.id = (' ', int(anrlist[i][1].strip()), ' ')
        i += 1
        
w = PDBIO()
w.set_structure(inpstruct)
w.save('1.pdb')
