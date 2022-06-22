#!/usr/bin/env python

from sys import argv
import os
import urllib
from Bio import SeqIO

script, fasta_file = argv

structure_dir="structures"

cwd = os.getcwd()

for record in SeqIO.parse(fasta_file, "fasta"):
    pdbid = record.id[:4]
    pdbfile = os.path.join(cwd,structure_dir,pdbid+".pdb.gz")
    if not os.path.isfile(pdbfile): 
        print pdbid, pdbfile
        urllib.urlretrieve ("https://files.rcsb.org/download/"+pdbid+".pdb.gz", "structures/"+pdbid+".pdb.gz")
