#! /usr/bin/env python3
import os
from Bio import SeqIO
import re
import json
from sys import argv
script, filename = argv

genemapdata = {}
subgroup_name, gene_name, allele_name, trseq, fullseq, halfseq, species = [], [], [], [], [], [], []
handle = open(filename, "rU")
for record in SeqIO.parse(handle, "fasta") :
    inseq = str(record.seq)
    if '*' in inseq:
        print inseq
        print "--"
        print record
        continue
    #split fasta header
    fastaheadersplit =  record.description.split("|")
    curr_allele_name = fastaheadersplit[1]
    curr_gene_name = curr_allele_name.split('*')[0]
    curr_subgroup_name = curr_gene_name.replace('/', '-').split('-')[0]
    curr_species = fastaheadersplit[2]
    allele_name.append(curr_allele_name)
    gene_name.append(curr_gene_name)
    subgroup_name.append(curr_subgroup_name)
    fullseq.append(str(inseq))
    species.append(str(curr_species))

genemapdata = [{"subgroup_name": a, "gene_name": b, "allele_name": c, "fullseq": d, "species":e} for a, b, c, d, e in zip(subgroup_name, gene_name, allele_name, fullseq, species)]

jsonfilename = os.path.splitext(filename)[0] + ".json"
print "writing..." ,jsonfilename
with open(jsonfilename, 'w') as fp:
    js = json.dumps(genemapdata)
    fp.write(js+"\n")
print "Done!"
