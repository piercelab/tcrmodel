#!/usr/bin/env python

import sys, os, re
import gzip
from sys import argv
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
from Bio import SeqIO

parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()

from TCR_functions import *
'''
script, pdbcode, chainida, chainidb = argv


gzpdbfile_path =  "/TCRmodeller/PDB_RELEASE/pdb_structures" +  '/%s/pdb%s.ent.gz' %(pdbcode[1:3], pdbcode) 
#gzpdbfile_path =  pdbcode+".pdb.gz" 
#import wget
#url = "http://www.rcsb.org/pdb/files/"+gzpdbfile_path
#filename = wget.download(url)
gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
pdbfile = parser.get_structure("PDB", gzpdbfile)

mychaina = pdbfile[0][chainida]
io.set_structure(mychaina)
io.save('tmpa.pdb', nonHetSelect())
achainseq = ""
for ppe in ppb.build_peptides(mychaina):
    achainseq += str(ppe.get_sequence())

mychainb = pdbfile[0][chainidb]
io.set_structure(mychainb)
io.save('tmpb.pdb', nonHetSelect())
bchainseq = ""
for ppe in ppb.build_peptides(mychainb):
    bchainseq += str(ppe.get_sequence())

fasta_file =  "/TCRmodeller/PDB_RELEASE/pdb_sequences/pdb_seqres.txt"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    if record.id[:4] == pdbcode and record.id[5:6] == chainida:
        afasta = record.seq 
    if record.id[:4] == pdbcode and record.id[5:6] == chainidb:
        bfasta = record.seq

#print achainseq
#print afasta

#print bchainseq
#print bfasta


print "\ntcrdomain"


#tmp1a = find_variable_domain_using_regex(achainseq, 'a')
tmp1a = assign_CDRs_using_REGEX(achainseq, 'A')
tmp1agr = tmp1a.groups()
#tmp2a = find_variable_domain_using_regex(afasta, 'a')
tmp2a = assign_CDRs_using_REGEX(afasta, 'A')
tmp2agr = tmp2a.groups()
 
#tmp1b = find_variable_domain_using_regex(bchainseq, 'b')
tmp1b = assign_CDRs_using_REGEX(bchainseq, 'B')
tmp1bgr = tmp1b.groups()
tmp2b = assign_CDRs_using_REGEX(bchainseq, 'B')
#tmp2b = find_variable_domain_using_regex(bfasta, 'b')
tmp2bgr = tmp2b.groups()

if re.search(tmp1agr[0], tmp1agr[0]) and re.search(tmp1bgr[0], tmp2bgr[0]):
    print "RES MATCH", pdbcode
else:
    print "RES NOMATCH", pdbcode
'''

fasta_file =  "tmp.fa"
#fasta_file =  "tcr_alpha.fasta"
#fasta_file =  "/TCRmodeller/PDB_RELEASE/pdb_sequences/pdb_seqres.txt" 
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    tmp1a = assign_CDRs_using_REGEX(record.seq, 'A')
    print tmp1a.groups()
#    if tmp1a is None:
    #print record.id, "FAIL"
