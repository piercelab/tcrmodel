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


import sys
sys.path.insert(0,'/www/tcrmodel/update/tcr_template_structures/scripts')

from TCR_functions import *
script, infile = argv


def cap_chain(pdbfile,capbegin,capend,chainid=None,outfile=None):
    st = parser.get_structure('PDB', pdbfile)
    i = 0#pdb_counter
    for chain in st[0]:
        if chainid is None:
            pass
        elif chain.id != chainid:
            continue
        for residue in list(chain):
            if (i < int(capbegin)) or (i >= int(capend)):
                chain.detach_child(residue.id)
            i += 1
    io.set_structure(st)
    if outfile is None:
        outfile = pdbfile+"_cap.pdb"
    io.save(outfile)
    return

cap_begin = 20
cap_end = 10

inpdb = infile+ ".trunc.fit.pdb"
#print inpdb

TAG = 'A'
chainchar = "D"
tmpcap_pdbfile = "tmpcap_"+inpdb
cap_pdbfile = infile + "_tcrpmhc.pdb"

fseq = pdbchain_to_fasta(inpdb, chainchar)
anarci_fseq = None 
if check_tcr_domain(fseq) is True:
    anarci_fseq = get_tcr_domain_seq(fseq, None)
    #print "here0", anarci_fseq
else:
    print inpdb, chainchar, "WARNING: tcr_alpha_beta_domain not found with ANARCI"
if anarci_fseq:
    regexres = assign_CDRs_using_REGEX(fseq, TAG)
    if regexres:
        #print regexres.groups()
        seq_vdomain = regexres.groups()[0]
        a_seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        pdbregex = re.search(seq_vdomain,fseq)
        if pdbregex:
            capbegin = pdbregex.span()[0]
            capend = pdbregex.span()[1]
            #print capbegin, capend, seq_vdomain, fseq
            cap_chain(inpdb, capbegin, capend, chainchar, tmpcap_pdbfile)
            
TAG = 'B'
chainchar = "E"
fseq = pdbchain_to_fasta(inpdb, chainchar)
anarci_fseq = None
if check_tcr_domain(fseq) is True:
    anarci_fseq = get_tcr_domain_seq(fseq, None)
    #print "here1", anarci_fseq
else:
    print pdbid, chainchar, "WARNING: tcr_alpha_beta_domain not found with ANARCI"
if anarci_fseq:
    regexres = assign_CDRs_using_REGEX(fseq, TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        b_seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        pdbregex = re.search(seq_vdomain,fseq)
        if pdbregex:
            capbegin = pdbregex.span()[0]
            capend = pdbregex.span()[1]
            cap_chain(tmpcap_pdbfile, capbegin, capend, chainchar, cap_pdbfile)


if os.path.isfile(cap_pdbfile): 
    pseq = pdbchain_to_fasta(inpdb, "C")
    plen= len(pseq)

    pf1 = open("MHC2a.fasta", "a+")
    pf1.write(">"+infile+"\n")
    mhc2aseq = pdbchain_to_fasta(inpdb, "A")
    pf1.write(mhc2aseq+"\n")

    pf2 = open("MHC2b.fasta", "a+")
    pf2.write(">"+infile+"\n")
    mhc2bseq = pdbchain_to_fasta(inpdb, "B")
    pf2.write(mhc2bseq+"\n")

    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+a_seq_fw+" "+b_seq_fw
    #print infile+" -mhc1 "+pdbchain_to_fasta(inpdb, "A")+" - peptide "+pdbchain_to_fasta(inpdb, "C")+" -alpha "+pdbchain_to_fasta(inpdb, "D")+" -beta "+pdbchain_to_fasta(inpdb, "E") + " -include_list in_list.txt -ignore_list ig_list.txt"
    print infile+"\t"+pdbchain_to_fasta(inpdb, "A")+"\t"+pdbchain_to_fasta(inpdb, "B")+"\t"+pdbchain_to_fasta(inpdb, "C")+"\t"+pdbchain_to_fasta(inpdb, "D")+"\t"+pdbchain_to_fasta(inpdb, "E")


