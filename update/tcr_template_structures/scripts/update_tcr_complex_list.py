#!/usr/bin/env python                

                                                                                                                         
#CREATE ORIENTATION template database
#create fasta seq file and complex structures
#need the previously created fw template seq and renumbered aho files

#TCR_ORIENTATION.seq
import os

import sys
sys.path.insert(0,'/home/rg/scripts')

from Bio import SeqIO
import glob
from TCR_functions import *



def seq_from_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id == inp_id:
            return rec.seq

outfile="TCR_complex.seq"
if os.path.exists(outfile):
    os.remove(outfile)


#first map ids of beta fw seqs
b_fw_fasta_file = "B_TCR_VDOMAIN_cap.fs"
b_fw = SeqIO.parse(open(b_fw_fasta_file),'fasta')
b_ids = map(lambda x: x.id, b_fw)

#loop thru alpha fw seqs and check matching ids
a_fw_fasta_file = "A_TCR_VDOMAIN_cap.fs"
a_fw = SeqIO.parse(open(a_fw_fasta_file),'fasta')
ori_list = []
for arecord in a_fw:
    apdb = arecord.id[:8]+"_aho.pdb"
    for brecord in b_ids:
        if (arecord.id[:4] == brecord[:4]):
            bpdb = brecord[:8]+"_aho.pdb"
            print apdb,bpdb
            dist = calc_alphacys_betacys_distance(apdb,bpdb)
            if (dist < 25):
                print apdb, bpdb, dist
                #complexpdb = make_complex(apdb,bpdb)
                #print complexpdb
                apdb_fw_seq = seq_from_id(arecord.id, a_fw_fasta_file)
                bpdb_fw_seq = seq_from_id(brecord, b_fw_fasta_file)
                ori_list.append([arecord.id[:8],brecord[:8],str(apdb_fw_seq),str(bpdb_fw_seq),str(brecord[9:])])
                with open(outfile, 'a') as fr:
                    fr.write(arecord.id[:8]+" "+brecord[:8]+" "+str(apdb_fw_seq)+" "+str(bpdb_fw_seq)+" "+str(brecord[9:])+"\n")



print "Hi"
#print ori_list
tmp_ori_list = ori_list
unq_ori_list = []
for item1 in ori_list:
    flag=True
    for unqitem in unq_ori_list:
        print "unqitem", unqitem
        print item1
        if ( str(item1[2])+str(item1[3]) == str(unqitem[2])+str(unqitem[3]) ):
            flag=False
    if flag:
        newitem = item1
        for item2 in tmp_ori_list:
            if ( str(newitem[2])+str(newitem[3]) == str(item2[2])+str(item2[3]) ):
                if (newitem[4] > item2[4]):
                    newitem = item2
        unq_ori_list.append(newitem)

uniqfile = "uniq.seq"
if os.path.exists(uniqfile):
    os.remove(uniqfile)

for col in unq_ori_list:
    with open(uniqfile, 'a') as fr:
        fr.write(col[0]+" "+col[1]+" "+col[2]+" "+col[3]+" "+col[4]+"\n")
        #make complex
        complexpdb = make_complex(col[0]+"_aho.pdb",col[1]+"_aho.pdb")

    TAG="A"
    outheader=str(col[0])+"_"+str(col[4])
    regexres = assign_CDRs_using_REGEX(col[2], TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        seq_gm = regexres.groups()[1] + regexres.groups()[2] + regexres.groups()[3] + regexres.groups()[4] + regexres.groups()[5] + regexres.groups()[7]
        seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seq_cdr1 = regexres.groups()[2]
        seq_cdr2hv4 = regexres.groups()[4]
        seq_cdr3 = regexres.groups()[6]
        with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_vdomain+"\n")
        with open(TAG+'_TCR_GM.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_gm+"\n")
        with open(TAG+'_TCR_FW.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_fw+"\n")
        with open(TAG+'_TCR_CDR1_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr1+"\n")
        with open(TAG+'_TCR_CDR2HV4_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr2hv4+"\n")
        with open(TAG+'_TCR_CDR3_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr3+"\n")

    TAG="B"
    outheader=str(col[1])+"_"+str(col[4])
    regexres = assign_CDRs_using_REGEX(col[3], TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        seq_gm = regexres.groups()[1] + regexres.groups()[2] + regexres.groups()[3] + regexres.groups()[4] + regexres.groups()[5] + regexres.groups()[7]
        seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seq_cdr1 = regexres.groups()[2]
        seq_cdr2hv4 = regexres.groups()[4]
        seq_cdr3 = regexres.groups()[6]
        with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_vdomain+"\n")
        with open(TAG+'_TCR_GM.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_gm+"\n")
        with open(TAG+'_TCR_FW.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_fw+"\n")
        with open(TAG+'_TCR_CDR1_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr1+"\n")
        with open(TAG+'_TCR_CDR2HV4_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr2hv4+"\n")
        with open(TAG+'_TCR_CDR3_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr3+"\n")

