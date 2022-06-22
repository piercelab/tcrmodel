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
#sys.path.insert(0,'~/tcr_update/scripts')
from TCR_functions import *
#script, inpdb = argv


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

cap_begin = 20
cap_end = 10

tcrpmhc_class1 = ["1ao7","1bd2","1fo0","1g6r","1kj2","1lp9","1mi5","1mwa","1nam","1oga","1qrn","1qse","1qsf","2ak4","2bnq","2bnr","2e7l","2esv","2f53","2f54","2gj6","2j8u","2jcc","2nx5","2oi9","2ol3","2p5e","2p5w","2pye","2uwe","2vlj","2vlk","2vlr","2ypl","3d39","3d3v","3dxa","3e2h","3e3q","3ffc","3gsn","3h9s","3hg1","3kpr","3kps","3kxf","3mv7","3mv8","3mv9","3o4l","3pqy","3pwp","3qdg","3qdj","3qdm","3qeq","3qfj","3rgv","3sjv","3tf7","3tfk","3tjh","3tpu","3uts","3utt","3vxm","3vxr","3vxs","3vxu","4eup","4ftv","4g8g","4g9f","4jfd","4jfe","4jff","4jrx","4jry","4l3e","4mji","4mnq","4ms8","4mvb","4mxq","4n0c","4n5e","4nhu","4prh","4pri","4prp","4qok","4qrp","4qrr","5brz","5bs0","5c07","5c08","5c09","5c0a","5c0b","5c0c","5d2l","5d2n","5e6i","5e9d","5eu6","5euo","5hhm","5hho","5hyj","5isz","5ivx","5jhd","5jzi","5m00","5m01","5m02","5men","5nht","5nme","5nmf","5nmg","5nqk","5sws","5swz","5tez","5til","5tje","5w1v","5w1w","5wkf","5wkh","5wlg","6am5","6amu","6avf","6avg","6bj2","6bj3","6bj8","6eqb","6mtm","6eqa","5yxu","5yxn","6g9q","6dkp","6d78","5xov","5xot","2ckb","6q3s"]

tcrpmhc_class2 = ["1D9K","1FYT","1J8H","1U3H","1YMM","1ZGL","2IAM","2IAN","2PXY","2WBJ","2Z31","3C5Z","3C60","3C6L","3MBE","3O6F","3PL6","3QIB","3QIU","3QIW","3RDT","4E41","4GG6","4GRL","4H1L","4MAY","4OZF","4OZG","4OZH","4OZI","4P23","4P2O","4P2Q","4P2R","4P46","4P4K","4P5T","4Y19","4Y1A","4Z7U","4Z7V","4Z7W","5KS9","5KSA","5KSB","6BGA","6CQL","6CQN","6CQQ","6CQR","6DFX","6DFW","6DFS","3T0E"]

fasta_file = "/www/tcrmodel/update/tcr_template_structures/rcsb_update/tcr_beta.fasta"
fhandle = open(fasta_file, "rU")
flist = list(SeqIO.parse(fhandle, "fasta"))
TAG = 'B'
chainchar = "E"
#tmpcap_pdbfile = "tmpcap_"+inpdb
pdbfile = "tmp.pdb"

for class1pdb in tcrpmhc_class2:
    for val in flist:
        pdbid = val.id[:4]
        if ((class1pdb.lower()) == pdbid):
            chainid = val.id[5:6]
            resolution = extract_resolution_from_pdb(val.id)
            pdbgzfile = os.path.join("/www/tcrmodel/update/tcr_template_structures/rcsb_update/structures/"+pdbid+".pdb.gz")
            get_pdb_from_gzpdbfile(val.id[:4],val.id[5:6],pdbfile,pdbgzfile)
            fseq = val.seq
            fseq = str(fseq).replace('X', '')
            pseq = pdb_to_fasta(pdbfile)
            #print pdbid, chainid
            #print fseq
            #print pseq
            fseq_vdomain = None
            pseq_vdomain = None


            fregexres = assign_CDRs_using_REGEX(fseq, TAG)
            if fregexres:
                fseq_vdomain = fregexres.groups()[0]
                seq_cdr1 = fregexres.groups()[2]
                seq_cdr2hv4 = fregexres.groups()[4]
                seq_cdr3 = fregexres.groups()[6]   
                #print pdbid, fseq_vdomain, seq_cdr1, seq_cdr2hv4, seq_cdr3
            else:
                pass
                #print pdbid, "fseq fail"

            pregexres = assign_CDRs_using_REGEX(pseq, TAG)
            if pregexres:
                pseq_vdomain = pregexres.groups()[0]
                seq_cdr1 = pregexres.groups()[2]
                seq_cdr2hv4 = pregexres.groups()[4]
                seq_cdr3 = pregexres.groups()[6]   
                #print class1pdb, pseq_vdomain, seq_cdr1, seq_cdr2hv4, seq_cdr3
            else:
                pass
                #print pdbid, "pseq fail"

                #pregexres = assign_CDRs_using_REGEX(pseq, TAG, True)
                #pseq_vdomain = pregexres.groups()[0]

            cdrfromseq = get_cdr_from_seq_by_aho_num_ext(fseq,TAG)
            seq_cdr1_aho = str(cdrfromseq[0])
            seq_cdr2_aho = str(cdrfromseq[1])  
            seq_cdr3_aho = str(cdrfromseq[4])  
            
            if ( (fregexres and pregexres) and (fseq_vdomain == pseq_vdomain) ):
                #print class1pdb, chainid, resolution, seq_cdr1_aho, seq_cdr2_aho, seq_cdr3_aho
                print class1pdb, resolution, seq_cdr1_aho, seq_cdr2_aho, seq_cdr3_aho
            else:
                #print class1pdb, chainid, resolution, seq_cdr1_aho, seq_cdr2_aho, seq_cdr3_aho, "no match"
                print class1pdb, resolution, seq_cdr1_aho, seq_cdr2_aho, seq_cdr3_aho, "no match"
            #break

