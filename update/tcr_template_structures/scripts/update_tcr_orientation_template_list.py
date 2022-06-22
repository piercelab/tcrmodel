#!/usr/bin/env python                

                                                                                                                         
#CREATE ORIENTATION template database
#create fasta seq file and complex structures
#need the previously created fw template seq and renumbered aho files

#TCR_ORIENTATION.seq
import os

import sys
sys.path.insert(0,'/www/tcrmodel/update/tcr_template_structures/scripts')

from Bio import SeqIO
import glob
from TCR_functions import *



def seq_from_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id == inp_id:
            return rec.seq

outfile="TCR_VD_ORIENTATION.seq"
if os.path.exists(outfile):
    os.remove(outfile)
no_ori_file="TCR_VD_NO_ORIENTATION.txt"
if os.path.exists(no_ori_file):
    os.remove(no_ori_file)

#first map ids of beta fw seqs
b_fw_fasta_file = "B_TCR_VDOMAIN.fasta"
b_fw = SeqIO.parse(open(b_fw_fasta_file),'fasta')
b_ids = map(lambda x: x.id, b_fw)

#loop thru alpha fw seqs and check matching ids
a_fw_fasta_file = "A_TCR_VDOMAIN.fasta"
a_fw = SeqIO.parse(open(a_fw_fasta_file),'fasta')
ori_list = []
for arecord in a_fw:
    apdb = arecord.id[:8]+"_aho.pdb"
    print apdb
    ori_match = False
    for brecord in b_ids:
        if (arecord.id[:4] == brecord[:4]):
            bpdb = brecord[:8]+"_aho.pdb"
            dist = calc_alphacys_betacys_distance(apdb,bpdb)
            if (dist < 25):
                ori_match = True
                print apdb, bpdb, dist
                complexpdb = make_complex(apdb,bpdb)
                #print complexpdb
                apdb_fw_seq = seq_from_id(arecord.id, a_fw_fasta_file)
                bpdb_fw_seq = seq_from_id(brecord, b_fw_fasta_file)
                ori_list.append([arecord.id[:8],brecord[:8],str(apdb_fw_seq),str(bpdb_fw_seq),str(brecord[9:])])
                with open(outfile, 'a') as fr:
                    fr.write(arecord.id[:8]+" "+brecord[:8]+" "+str(apdb_fw_seq)+" "+str(bpdb_fw_seq)+" "+str(brecord[9:])+"\n")
    if not ori_match:
        with open(no_ori_file, 'a') as fr:
            fr.write(arecord.id[:8]+"\n")
        
