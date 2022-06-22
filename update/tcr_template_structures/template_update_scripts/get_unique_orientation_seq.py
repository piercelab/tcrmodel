#!/usr/bin/env python

import os
from glob import glob
from Bio import SeqIO
from update_functions import *
import shutil

ori_seq_file="TCR_VD_ORIENTATION.seq"
with open(ori_seq_file, 'r') as f:
    ori_list = [ line.split() for line in f ]
#
#Find unique Non-redundant orientation seq    
#
nr_ori_seq_file = "NR_TCR_VD_ORIENTATION.seq"
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
#
#Find redundant seqs list in orientation 
#to remove from other segments fasta files
#
AB_ori_redundant_file = "TCR_VD_ORIENTATION_Redundant_list.txt"
differences = []
with open(AB_ori_redundant_file, 'w') as rd:
    for fitem in ori_list:
        if fitem not in unq_ori_list:
            differences.append(fitem)
            rd.write(fitem[0]+"\n")
            rd.write(fitem[1]+"\n")


#
#Find redundant seqs in sequences that are not oriented(no complex formed)
#to remove from other segments fasta files
#
with open(nr_ori_seq_file, 'r') as f:
    nr_ori_list = [ line.split() for line in f ]

A_seq_clean_file = "A_TCR_VDOMAIN_clean.fasta"
A_id_file = "A_TCR_VD_NO_ORIENTATION.txt"
with open(A_id_file, 'r') as f:
    A_id_list = [ line.strip() for line in f ]
A_no_ori_redundant_file_cmplx = "A_TCR_VD_NO_ORIENTATION_R_cmplx.txt"
with open(A_no_ori_redundant_file_cmplx, 'w') as af:
    for myid in A_id_list:
        myseq = seq_from_id(myid,A_seq_clean_file)
        for col in nr_ori_list:
            if ( str(col[2].strip()) == str(myseq.strip()) ):
                af.write(myid+"\n")
                break

B_seq_clean_file = "B_TCR_VDOMAIN_clean.fasta"
B_id_file = "B_TCR_VD_NO_ORIENTATION.txt"
with open(B_id_file, 'r') as f:
    B_id_list = [ line.strip() for line in f ]
B_no_ori_redundant_file_cmplx = "B_TCR_VD_NO_ORIENTATION_R_cmplx.txt"
with open(B_no_ori_redundant_file_cmplx, 'w') as bf:
    for myid in B_id_list:
        myseq = seq_from_id(myid,B_seq_clean_file)
        for col in nr_ori_list:
            if ( str(col[3].strip()) == str(myseq.strip()) ):
                bf.write(myid+"\n")
                break

A_no_ori_redundant_file_self = "A_TCR_VD_NO_ORIENTATION_R_self.txt"
B_no_ori_redundant_file_self = "B_TCR_VD_NO_ORIENTATION_R_self.txt"
unwanted_list_of_files = [
    AB_ori_redundant_file, 
    A_no_ori_redundant_file_cmplx,
    B_no_ori_redundant_file_cmplx,
    A_no_ori_redundant_file_self,
    B_no_ori_redundant_file_self
]

unwanted_seqs = set()
for unwanted_file in unwanted_list_of_files:
    with open(unwanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                unwanted_seqs.add(line[:8])



extensions = (
    'A_TCR_GM.fasta', 
    'B_TCR_GM.fasta',
    'A_TCR_FW.fasta', 
    'B_TCR_FW.fasta',
    'A_TCR_CDR*.fasta',
    'B_TCR_CDR*.fasta'
)
curr_list = []
for extension in extensions:
    curr_list.extend(glob(extension))
for curr_file in curr_list:
    fasta_sequences = SeqIO.parse(open(curr_file),'fasta')
    result_file = "NR1_"+ curr_file
    with open(result_file, "w") as f:
        for seq in fasta_sequences:
            if seq.id[:8] in unwanted_seqs:
                continue
            SeqIO.write([seq], f, "fasta")


for curr_file in curr_list:
    nr_file = "NR1_"+ curr_file
    unq_list = get_unique_by_seq_and_resolution(nr_file)
    nr_file = "NR2_"+ curr_file
    if unq_list:
        with open(nr_file, "w") as f:
            for dat in unq_list:
                SeqIO.write(dat[3], f, "fasta")


#
#remove redundancy from other sequences
#
#curr_list=['B_TCR_CDR2_11.fasta']
for curr_file in curr_list:
    other_file = "other_"+ curr_file
    nr_file = "NR2_"+ curr_file
    nr_other_file = "NR1_"+ other_file
    join_file =  "NR_"+ curr_file

    if os.path.isfile(other_file):
        print curr_file, other_file
        unq_list = get_unique_by_seq_and_resolution(other_file)
        if unq_list:
            with open(nr_other_file, "w") as f:
                for dat in unq_list:
                    fasta_file = nr_file 
                    seq_match = check_seq_match(dat[1],fasta_file)
                    if not seq_match:
                        SeqIO.write(dat[3], f, "fasta")
            with open(join_file,'wb') as wfd:
                for f in [nr_file, nr_other_file]:
                    with open(f,'rb') as fd:
                        shutil.copyfileobj(fd, wfd)

    else:
        shutil.copyfile(nr_file, join_file)


#
#combine all NR pdb chains from single PBD structure
#

pdb_achain_list = []
for curr_file in glob("NR_A_*.fasta"):
    fasta_sequences = SeqIO.parse(open(curr_file),'fasta')
    for val in fasta_sequences:
        pdb_achain_list.append(val.id[:8])
for curr_file in glob("NR_*.seq"):
    with open(curr_file,'r') as r:
        line_list = [ line.split() for line in r ]
        for line in line_list:
            pdb_achain_list.append(line[0])
list_a_chain = list(set(pdb_achain_list))


pdb_bchain_list = []
for curr_file in glob("NR_B_*.fasta"):
    fasta_sequences = SeqIO.parse(open(curr_file),'fasta')
    for val in fasta_sequences:
        pdb_bchain_list.append(val.id[:8])
for curr_file in glob("NR_*.seq"):
    with open(curr_file,'r') as r:
        line_list = [ line.split() for line in r ]
        for line in line_list:
            pdb_bchain_list.append(line[1])
list_b_chain = list(set(pdb_bchain_list))


a_nr_vd_file = "NR_A_TCR_VD.fasta"
a_vd_all_file = "A_TCR_VD_ALL.fasta"
with open(a_nr_vd_file, 'w') as rd:
    for chain in list_a_chain:
        print chain, "hi"
        curr_rec = fasta_record_from_trunc_id(chain, a_vd_all_file)
        if curr_rec:
            SeqIO.write(curr_rec, rd, "fasta")

b_nr_vd_file = "NR_B_TCR_VD.fasta"
b_vd_all_file = "B_TCR_VD_ALL.fasta"
with open(b_nr_vd_file, 'w') as rd:
    for chain in list_b_chain:
        curr_rec = fasta_record_from_trunc_id(chain, b_vd_all_file)
        print chain, curr_rec
        SeqIO.write(curr_rec, rd, "fasta")

myunqlist = list_a_chain + list_b_chain
counted = []
pdbchains_list = []
mytmplist = myunqlist
for val in myunqlist:
    if not val[:4] in counted:
        pdbchains = []
        for tmpval in mytmplist:
            if val[:4] == tmpval[:4]:
                pdbchains.append(tmpval)
        counted.append(val[:4])
        pdbchains_list.append(pdbchains)

for pdbchains in pdbchains_list:
    make_complex_from_chains(pdbchains)

for pdbchains in pdbchains_list:
    print pdbchains
    if len(pdbchains) < 2:
        continue
    else:
        atag = pdbchains[0]
        btag = pdbchains[1]
        achainid = atag[5:6]
        bchainid = btag[5:6]
        if (achainid == bchainid):
            if achainid.isupper(): newbchainid = bchainid.lower()
            if achainid.islower(): newbchainid = bchainid.upper()
            new_btag = btag[:5] + newbchainid + btag[6:]
            for curr_file in glob("NR_*.fasta"):
                with open(curr_file, 'r') as f:
                    filedata = f.read()
                filedata = filedata.replace(btag, new_btag)
                with open(curr_file, 'w') as f:
                    f.write(filedata+'\n')
                    #f.write('\n') 

for curr_file in glob("NR_*.fasta"):
    cwd = os.getcwd()
    seq_dirpath = os.path.join(cwd,"../updated_templates/tcr/seq")
    if not os.path.exists(seq_dirpath):
        os.makedirs(seq_dirpath)
    new_file = os.path.join(seq_dirpath,curr_file[3:])
    shutil.copyfile(curr_file, new_file)

for curr_file in glob("NR_*.seq"):
    cwd = os.getcwd()
    seq_dirpath = os.path.join(cwd,"../updated_templates/tcr/seq")
    if not os.path.exists(seq_dirpath):
        os.makedirs(seq_dirpath)
    new_file = os.path.join(seq_dirpath,curr_file[3:])
    shutil.copyfile(curr_file, new_file)


