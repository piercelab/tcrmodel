#!/usr/bin/env python

import os
ori_seq_file="TCR_VD_ORIENTATION.seq"

with open(ori_seq_file, 'r') as f:
    ori_list = [ line.split() for line in f ]
    
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

uniqfile = "TCR_VD_ORIENTATION_NR.seq"
if os.path.exists(uniqfile):
    os.remove(uniqfile)
with open(uniqfile, 'a') as fr:
    for col in unq_ori_list:
        fr.write(col[0]+" "+col[1]+" "+col[2]+" "+col[3]+" "+col[4]+"\n")


redundantfile = "TCR_VD_ORIENTATION_Redundant_list.txt"
differences = []

with open(redundantfile, 'a') as rd:
    for fitem in ori_list:
        if fitem not in unq_ori_list:
            differences.append(fitem)
            rd.write(fitem[0]+"\n")
            rd.write(fitem[1]+"\n")


print differences
