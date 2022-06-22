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
# sys.path.insert(0,'~/tcr_update/scripts')
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

def fasta_record_from_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id[:4] == inp_id:
            return rec.seq

#print infile
cap_begin = 20
cap_end = 10
TAG = 'A'
chainchar = "D"
#inpdb = infile+ ".trunc.fit.pdb"
inpdb = infile + "_tcrpmhc.pdb"
tmpcap_pdbfile = "tmpcap_"+inpdb
cap_pdbfile = infile + "_tcrpmhc.pdb"
a_fasta_file = "/www/tcrmodel/update/tcr_template_structures/rcsb_update/tcr_alpha.fasta"
d_fasta_file = "/www/tcrmodel/update/tcr_template_structures/rcsb_update/tcr_delta.fasta"
a_fasta_seq = fasta_record_from_id(infile, a_fasta_file)
if not a_fasta_seq:
    a_fasta_seq = fasta_record_from_id(infile, d_fasta_file)
a_fasta_seq = str(a_fasta_seq).replace('X', '') 

b_fasta_file = "/www/tcrmodel/update/tcr_template_structures/rcsb_update/tcr_beta.fasta"
g_fasta_file = "/www/tcrmodel/update/tcr_template_structures/rcsb_update/tcr_gamma.fasta"
b_fasta_seq = fasta_record_from_id(infile, b_fasta_file)
if not b_fasta_seq:
    b_fasta_seq = fasta_record_from_id(infile, g_fasta_file)
b_fasta_seq = str(b_fasta_seq).replace('X', '') 

a_pdb_seq = pdbchain_to_fasta(inpdb, "D")
b_pdb_seq = pdbchain_to_fasta(inpdb, "E")

#print "achain"
#print a_fasta_seq
#print a_pdb_seq
#print get_capped_tcr_vdomain(a_fasta_seq)
#print get_capped_tcr_vdomain(a_pdb_seq)
#print "bchain"
#print b_fasta_seq
#print b_pdb_seq
#print get_capped_tcr_vdomain(b_fasta_seq)
#print get_capped_tcr_vdomain(b_pdb_seq)
a_fa_cap = get_capped_tcr_vdomain(a_fasta_seq)
a_pdb_cap = get_capped_tcr_vdomain(a_pdb_seq)
b_fa_cap = get_capped_tcr_vdomain(b_fasta_seq)
b_pdb_cap = get_capped_tcr_vdomain(b_pdb_seq)
'''
print a_fa_cap
print a_pdb_seq
print a_pdb_cap
print b_fa_cap
print b_pdb_seq
print b_pdb_cap
'''
info = get_pdb_info(infile)
print "HI:", infile, get_capped_tcr_vdomain(a_fasta_seq), get_capped_tcr_vdomain(b_fasta_seq), info[0], info[1] 
#if (a_pdb_seq in a_fasta_seq) and (b_pdb_seq in b_fasta_seq):
if (a_pdb_cap == a_fa_cap) and (b_pdb_cap == b_fa_cap):
    pass
    #print infile, "yes", info[0], info[1]
else:
    pass
    #print infile, "no", info[0], info[1]



'''
a_f_segs = get_cdr_from_seq_by_aho_num(get_capped_tcr_vdomain(a_fasta_seq), "A", True)
#a_f_segs = get_cdr_from_seq_by_aho_num(a_fasta_seq, "A", True)
a_p_segs = get_cdr_from_seq_by_aho_num(get_capped_tcr_vdomain(a_pdb_seq), "A", True)
#a_p_segs = get_cdr_from_seq_by_aho_num(a_pdb_seq, "A", True)
b_f_segs = get_cdr_from_seq_by_aho_num(get_capped_tcr_vdomain(b_fasta_seq), "B", True)
b_p_segs = get_cdr_from_seq_by_aho_num(get_capped_tcr_vdomain(b_pdb_seq), "B", True)

if ( (a_f_segs[5] == a_p_segs[5]) and (b_f_segs[5] == b_p_segs[5]) ):
    print infile, "YES"
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+a_f_segs[5]+" "+b_f_segs[5]
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+a_p_segs[5]+" "+b_p_segs[5]
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+get_capped_tcr_vdomain(a_fasta_seq)+" "+get_capped_tcr_vdomain(b_fasta_seq)
    #cap_pdb_chain_segment(inpdb, "D", 4, 148, tmpcap_pdbfile)
    #cap_pdb_chain_segment(tmpcap_pdbfile, "E", 4, 148, cap_pdbfile)
    #print a_f_segs
    #print a_p_segs
    #print b_f_segs
    #print b_p_segs
    #pass
else:
    print infile, "NO"
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+a_p_segs[5]+" "+b_p_segs[5]
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "C")+" "+get_capped_tcr_vdomain(a_fasta_seq)+" "+get_capped_tcr_vdomain(b_fasta_seq)
    #print a_f_segs
    #print a_p_segs
    #print b_f_segs
    #print b_p_segs
'''
'''
TAG="A"
chainchar="D"
anarci_apseq = None 
if check_tcr_domain(a_pdb_seq) is True:
    anarci_apseq = get_tcr_domain_seq(a_pdb_seq, None)
else:
    print inpdb, chainchar, "WARNING: tcr_alpha_beta_domain not found with ANARCI"

if anarci_apseq:
    regexres = assign_CDRs_using_REGEX(a_pdb_seq, TAG)
    if regexres:
        #print regexres.groups()
        seq_vdomain = regexres.groups()[0]
        a_seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        pdbregex = re.search(seq_vdomain,a_pdb_seq)
        if pdbregex:
            capbegin = pdbregex.span()[0]
            capend = pdbregex.span()[1]
            cap_chain(inpdb, capbegin, capend, chainchar, tmpcap_pdbfile)
        else:
            print "pdbregex failed!"
    else:
        print "regex failed!"
            
TAG = 'B'
chainchar = "E"
anarci_bpseq = None
if check_tcr_domain(b_pdb_seq) is True:
    anarci_bpseq = get_tcr_domain_seq(b_pdb_seq, None)
else:
    print pdbid, chainchar, "WARNING: tcr_alpha_beta_domain not found with ANARCI"
if anarci_bpseq:
    regexres = assign_CDRs_using_REGEX(b_pdb_seq, TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        b_seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        pdbregex = re.search(seq_vdomain,b_pdb_seq)
        if pdbregex:
            capbegin = pdbregex.span()[0]
            capend = pdbregex.span()[1]
            cap_chain(tmpcap_pdbfile, capbegin, capend, chainchar, cap_pdbfile)
'''
if os.path.isfile(cap_pdbfile): 
    pseq = pdbchain_to_fasta(inpdb, "C")
    plen= len(pseq)
    pfname = "MHC1a.P"+str(plen)+".fasta"
    pf=open(pfname, "a+")
    pf.write(">"+infile+"\n")
    mhcseq = pdbchain_to_fasta(inpdb, "A")
    pf.write(mhcseq+"\n")
    #print infile+" "+pdbchain_to_fasta(inpdb, "A")+" "+pdbchain_to_fasta(inpdb, "B")+" "+pdbchain_to_fasta(inpdb, "C")+" "+a_seq_fw+" "+b_seq_fw
    print infile+"\t"+pdbchain_to_fasta(inpdb, "A")+"\t"+pdbchain_to_fasta(inpdb, "B")+"\t"+pdbchain_to_fasta(inpdb, "C")+"\t"+a_fa_cap+"\t"+b_fa_cap
    #print infile+" -mhc1 "+pdbchain_to_fasta(inpdb, "A")+" - peptide "+pdbchain_to_fasta(inpdb, "C")+" -alpha "+pdbchain_to_fasta(inpdb, "D")+" -beta "+pdbchain_to_fasta(inpdb, "E") + " -include_list in_list.txt -ignore_list ig_list.txt"
    #print infile+"\t"+pdbchain_to_fasta(inpdb, "A")+"\t"+pdbchain_to_fasta(inpdb, "C")+"\t"+pdbchain_to_fasta(inpdb, "D")+"\t"+pdbchain_to_fasta(inpdb, "E")
    #print infile+"\t"+pdbchain_to_fasta(inpdb, "A")+"\t"+pdbchain_to_fasta(inpdb, "B")+"\t"+pdbchain_to_fasta(inpdb, "C")+"\t"+pdbchain_to_fasta(inpdb, "D")+"\t"+pdbchain_to_fasta(inpdb, "E")

