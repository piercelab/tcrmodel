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
from shutil import copyfile

from update_functions import *

outfile="TCR_VD_ORIENTATION.seq"
if os.path.exists(outfile):
    os.remove(outfile)

#first map ids of beta fw seqs
b_fw_fasta_file = "B_TCR_VDOMAIN_clean.fasta"
b_fw = SeqIO.parse(open(b_fw_fasta_file),'fasta')
b_ids = map(lambda x: x.id, b_fw)

#loop thru alpha fw seqs and check matching ids
a_fw_fasta_file = "A_TCR_VDOMAIN_clean.fasta"
a_fw = SeqIO.parse(open(a_fw_fasta_file),'fasta')

a_ori_list = []
b_ori_list = []
ori_list = []
a_all_list = []

for arecord in a_fw:
    a_all_list.append(arecord.id)
    apdb = arecord.id[:8]+"_aho.pdb"
    ori_match = False
    for brecord in b_ids:
        if (arecord.id[:4] == brecord[:4]):
            bpdb = brecord[:8]+"_aho.pdb"
            dist = calc_alphacys_betacys_distance(apdb,bpdb)
            if (dist < 25):
                ori_match = True
                #print apdb, bpdb, dist
                #complex_pdbname = brecord[:4] + "_aho.pdb"
                #print complex_pdbname
                #complexpdb = make_complex(apdb,bpdb,complex_pdbname)
                complexpdb = make_complex(apdb,bpdb)
                #copyfile(complexpdb, "tcr/pdb/"+complexpdb)
                #print complexpdb
                apdb_fw_seq = seq_from_id(arecord.id, a_fw_fasta_file)
                bpdb_fw_seq = seq_from_id(brecord, b_fw_fasta_file)
                ori_list.append([arecord.id[:8],brecord[:8],str(apdb_fw_seq),str(bpdb_fw_seq),str(brecord[9:])])
                a_ori_list.append(arecord.id)
                b_ori_list.append(brecord)
                with open(outfile, 'a') as fr:
                    fr.write(arecord.id[:8]+" "+brecord[:8]+" "+str(apdb_fw_seq)+" "+str(bpdb_fw_seq)+" "+str(brecord[9:])+"\n")
                    releaseDate_val, resolution_val = get_pdb_info(arecord.id[:4])
                    print arecord.id[:4].upper(), str(apdb_fw_seq), str(bpdb_fw_seq), releaseDate_val, resolution_val

b_all_list = []
for brecord in b_ids:
    b_all_list.append(brecord)

a_no_ori_list = []
a_no_ori_file = "A_TCR_VD_NO_ORIENTATION.txt"
with open(a_no_ori_file, 'w') as rd:
    for item in a_all_list:
        if item not in a_ori_list:
            rd.write(item+"\n")
            a_no_ori_list.append(item)
a_no_ori_seq = "A_TCR_VD_NO_ORIENTATION.fasta"
with open(a_no_ori_seq, 'w') as ard:
    for item in a_all_list:
        if item not in a_ori_list:
            curr_rec = fasta_record_from_id(item, a_fw_fasta_file)
            SeqIO.write(curr_rec, ard, "fasta")
a_no_ori_list_nr = []
unq_list = get_unique_by_seq_and_resolution(a_no_ori_seq)
if unq_list:
    nr_file = "NR1_"+ a_no_ori_seq
    with open(nr_file, "w") as f:
        for dat in unq_list:
            SeqIO.write(dat[3], f, "fasta")
            a_no_ori_list_nr.append(dat[0])

b_no_ori_list = []
b_no_ori_file = "B_TCR_VD_NO_ORIENTATION.txt"
with open(b_no_ori_file, 'w') as rd:
    for item in b_all_list:
        if item not in b_ori_list:
            rd.write(item+"\n")
            b_no_ori_list.append(item)
b_no_ori_seq = "B_TCR_VD_NO_ORIENTATION.fasta"
with open(b_no_ori_seq, 'w') as brd:
    for item in b_all_list:
        if item not in b_ori_list:
            curr_rec = fasta_record_from_id(item, b_fw_fasta_file)
            SeqIO.write(curr_rec, brd, "fasta")
b_no_ori_list_nr = []
unq_list = get_unique_by_seq_and_resolution(b_no_ori_seq)
if unq_list:
    nr_file = "NR1_"+ b_no_ori_seq
    with open(nr_file, "w") as f:
        for dat in unq_list:
            SeqIO.write(dat[3], f, "fasta")
            b_no_ori_list_nr.append(dat[0])


A_no_ori_redundant_file_self = "A_TCR_VD_NO_ORIENTATION_R_self.txt"
with open(A_no_ori_redundant_file_self, 'w') as rd:
    for dat in list(set(a_no_ori_list) - set(a_no_ori_list_nr)):
        rd.write(dat+"\n")

B_no_ori_redundant_file_self = "B_TCR_VD_NO_ORIENTATION_R_self.txt"
with open(B_no_ori_redundant_file_self, 'w') as rd:
    for dat in list(set(b_no_ori_list) - set(b_no_ori_list_nr)):
        rd.write(dat+"\n")
