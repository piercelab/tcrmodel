#!/usr/bin/env python                

                                                                                                                         
#CREATE ORIENTATION template database
#create fasta seq file and complex structures
#need the previously created fw template seq and renumbered aho files

#TCR_ORIENTATION.seq
import os
import shutil
import sys
sys.path.insert(0,'/www/tcrmodel/update/tcr_template_structures/scripts')

from Bio import SeqIO
import glob
from TCR_functions import *
from shutil import copyfile

from update_functions import *

outfile="R_TCR_FW_ORIENTATION.seq"
if os.path.exists(outfile):
    os.remove(outfile)

#first map ids of beta fw seqs
b_fw_fasta_file = "B_TCR_FW.fasta"
b_fw = SeqIO.parse(open(b_fw_fasta_file),'fasta')
b_ids = map(lambda x: x.id, b_fw)

#loop thru alpha fw seqs and check matching ids
a_fw_fasta_file = "A_TCR_FW.fasta"
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
                print apdb, bpdb, dist
                #complex_pdbname = brecord[:4] + "_aho.pdb"
                #print complex_pdbname
                #complexpdb = make_complex(apdb,bpdb,complex_pdbname)
                complexpdb = make_complex(apdb,bpdb)
                #copyfile(complexpdb, "tcr/pdb/"+complexpdb)
                print complexpdb
                apdb_fw_seq = seq_from_id(arecord.id, a_fw_fasta_file)
                bpdb_fw_seq = seq_from_id(brecord, b_fw_fasta_file)
                ori_list.append([arecord.id[:8],brecord[:8],str(apdb_fw_seq),str(bpdb_fw_seq),str(brecord[9:])])
                a_ori_list.append(arecord.id)
                b_ori_list.append(brecord)
                with open(outfile, 'a') as fr:
                    fr.write(arecord.id[:8]+" "+brecord[:8]+" "+str(apdb_fw_seq)+" "+str(bpdb_fw_seq)+" "+str(brecord[9:])+"\n")




with open(outfile, 'r') as f:
    ori_list = [ line.split() for line in f ]
#                                                                                                                                             
#Find unique Non-redundant orientation seq                                                                                                    
#                                                                                                                                             
nr_ori_seq_file = "NR_TCR_FW_ORIENTATION.seq"
tmp_ori_list = ori_list
unq_ori_list = []
for item1 in ori_list:
    flag=True
    for unqitem in unq_ori_list:
        if ( str(item1[2])+str(item1[3]) == str(unqitem[2])+str(unqitem[3]) ):
            flag=False
    if flag:
        newitem = item1
        for item2 in tmp_ori_list:
            if ( str(newitem[2])+str(newitem[3]) == str(item2[2])+str(item2[3]) ):
                if (newitem[4] > item2[4]):
                    newitem = item2
        unq_ori_list.append(newitem)
with open(nr_ori_seq_file, 'w') as fr:
    for col in unq_ori_list:
        fr.write(col[0]+" "+col[1]+" "+col[2]+" "+col[3]+" "+col[4]+"\n")

cwd = os.getcwd()
seq_dirpath = os.path.join(cwd,"../updated_templates/tcr/seq")
if not os.path.exists(seq_dirpath):
    os.makedirs(seq_dirpath)
new_file = os.path.join(seq_dirpath,"TCR_FW_ORIENTATION.seq")
shutil.copyfile(nr_ori_seq_file,new_file)
