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

script, fasta_file, TAG = argv
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    pdbid = record.id[:4]

    if (pdbid == "4zez"):
        pdbid = "5jzi"
    chainid = record.id[5:6]
    outtag = pdbid+"_"+chainid+"_"+TAG
    pdbfile = outtag+".pdb"
    pdb_release_path = "/TCRmodeller/PDB_RELEASE/pdb_structures"
    gzpdbfile_path =  pdb_release_path + '/%s/pdb%s.ent.gz' %(pdbid[1:3], pdbid)
    if os.path.isfile(gzpdbfile_path):
        gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
        st = parser.get_structure(pdbid, gzpdbfile)
        deposition_date = st.header['deposition_date']
        resolution = st.header['resolution']
        print pdbid, deposition_date
        with open('tcr_depdate.txt','a') as fr:
            fr.write(pdbid+" "+deposition_date+"\n")
    else:
        print pdbid, "Does not exists\n"
