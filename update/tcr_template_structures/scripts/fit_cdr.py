#!/usr/bin/python 

import os
import glob
import subprocess
from subprocess import Popen, PIPE



def create_profit_infile(infile):

    f = open(infile,"w")             
    f.write("ATOMS N,CA,C,O\n")                             
    f.write("ZONE CLEAR\n")
    f.write("ZONE A*:A*\n")
    f.write("FIT\n")
    f.close()                                               
    
def get_rms(infile, reference, mobile):                                                     
    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    cmd = subprocess.Popen([profit_program, '-f', infile, '-h', reference, mobile], stdout=subprocess.PIPE)
    for line in cmd.stdout:
        if "RMS:" in line:
            rmsd = line.split(' ')[-1]
            return rmsd

def gen_pairwise_rmsd_matrix(listfiles,ofile):
        
    f = open(ofile,'w')
    listfiles_csv = ",".join(listfiles)
    f.write(listfiles_csv+"\n")

    for i in listfiles:
        all_rms = []  
        for j in listfiles:
            rms = get_rms(profit_infile,i,j)
            all_rms.append(rms.strip())
        all_rms_csv = ",".join(all_rms)
        f.write(all_rms_csv+"\n")

    f.close()

def gen_pairwise_rmsd_matrix_with_color(listfiles,ofile):
        
    rf = open(ofile+".R",'w')
    f = open(ofile,'w')
    listfiles_csv = ",".join(listfiles)
    f.write(listfiles_csv+"\n")
    all_colors = []

    for i in listfiles:
        
        if i.endswith('_L3_13.pdb'):all_colors.append('"red"')
        elif i.endswith('_A_cdr3_13.pdb'):all_colors.append('"green"')

        all_rms = []  
        for j in listfiles:
            rms = get_rms(profit_infile,i,j)
            all_rms.append(rms.strip())
            
        all_rms_csv = ",".join(all_rms)
        all_colors_csv = ",".join(all_colors)
        f.write(all_rms_csv+"\n")

    f.close()
    rf.write("cluscolors=c("+all_colors_csv+")"+"\n")
    rf.close()

# main
myPath= "/home/rg/tcr_db_update/tcr_template_library"

profit_infile = 'profit.in'
create_profit_infile(profit_infile)



cdrpdbtypes = ['A_cdr1', 'A_cdr2', 'A_cdr3', 'A_hv4', 'B_cdr1', 'B_cdr2', 'B_cdr3', 'B_hv4']
cdrlentypes = range(1, 26)
for f in cdrpdbtypes:
    for t in cdrlentypes:
        ftag = "*_*_"+str(f)+"_"+str(t)+".pdb"
        files_grabbed = []
        files_grabbed.extend(glob.glob(ftag))
        if len(files_grabbed) > 0:
            ofile = str(f)+"_"+str(t)+".rms"
            print ftag, ofile
            gen_pairwise_rmsd_matrix(files_grabbed,ofile)

   
         

'''
cdrlentypes = range(1, 26)
for f in cdrlentypes:   
    types = ('*_B_cdr3_'+str(f)+'.pdb', '*_H3_'+str(f)+'.pdb')
    files_grabbed = []
    for files in types:
        files_grabbed.extend(glob.glob(files))
    print "types" ,types, files, len(files_grabbed)
    if len(files_grabbed) > 0:
        ofile = 'CDR3_H3_'+str(f)+'.rms'
        gen_pairwise_rmsd_matrix(files_grabbed,ofile)
'''
