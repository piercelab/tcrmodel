#!/usr/bin/env python

import sys, os, re
from sys import argv
from Bio import SeqIO

import sys
sys.path.insert(0,'/www/tcrmodel/update/tcr_template_structures/scripts')
from TCR_functions import *

script, fasta_file, TAG = argv

pdbid_listoflists = []
pdbid_list = []
unique_fasta_file = TAG+'.unique.fasta'
nocap_fasta_file = TAG+'.nocap.fasta'
if os.path.exists(unique_fasta_file):
    os.remove(unique_fasta_file)
if os.path.exists(nocap_fasta_file):
    os.remove(nocap_fasta_file)
    

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    pdbid = record.id[:4]
    chainid = record.id[5:6]
    outtag = pdbid+"_"+chainid+"_"+TAG
    pdbfile = outtag+".pdb"
    print pdbid
    if pdbid in pdbid_list:
        print pdbid,"in pdblist"
        continue

    pdbgzfile = os.path.join("/www/tcrmodel/update/tcr_template_structures/rcsb_update/structures/"+pdbid+".pdb.gz")
    if not os.path.isfile(pdbgzfile):
        print "FILE not found", pdbgzfile
    get_pdb_from_gzpdbfile(record.id[:4],record.id[5:6],pdbfile,pdbgzfile)
    strseq = pdb_to_fasta(pdbfile)
    anarci_fseq = get_tcr_domain_seq(record.seq,TAG)
    print anarci_fseq
    regexres = assign_CDRs_using_REGEX(anarci_fseq, TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
    else:
        regexres = assign_CDRs_using_REGEX(record.seq, TAG)
        if regexres:
            seq_vdomain = regexres.groups()[0]
        else:
            with open(nocap_fasta_file, 'a') as fr:
                SeqIO.write(record, fr, "fasta")
            print "no match found", pdbid,chainid
            continue
    pdbregex = re.search(seq_vdomain,strseq)
    if pdbregex:
        pdbid_list.append(pdbid)
        with open(unique_fasta_file, 'a') as fr:
            SeqIO.write(record, fr, "fasta")
    else:
        with open(nocap_fasta_file, 'a') as fr:
            SeqIO.write(record, fr, "fasta")
            

if os.path.exists(nocap_fasta_file):
    nocap_sequences = SeqIO.parse(open(nocap_fasta_file),'fasta')
    for record in nocap_sequences:
        pdbid = record.id[:4]
        if pdbid in pdbid_list:
            continue
        else:
            print "adding", pdbid
            pdbid_list.append(pdbid)
            with open(unique_fasta_file, 'a') as fr:
                SeqIO.write(record, fr, "fasta")
