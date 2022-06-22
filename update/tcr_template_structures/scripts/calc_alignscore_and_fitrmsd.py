#!/usr/bin/python 

import os
import glob
import subprocess
from subprocess import Popen, PIPE

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocol -ignore_zero_occupancy false -ignore_unrecognized_res")


def calc_alignment_score(cdr_loop_seq, template_loop_seq):
    
    PAM30 = ([6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2],
             [-7,8,-6,-10,-8,-2,-9,-9,-2,-5,-8,0,-4,-9,-4,-3,-6,-2,-10,-8],
             [-4,-6,8,2,-11,-3,-2,-3,0,-5,-7,-1,-9,-9,-6,0,-2,-8,-4,-8],
             [-3,-10,2,8,-14,-2,2,-3,-4,-7,-12,-4,-11,-15,-8,-4,-5,-15,-11,-8],
             [-6,-8,-11,-14,10,-14,-14,-9,-7,-6,-15,-14,-13,-13,-8,-3,-8,-15,-4,-6],
             [-4,-2,-3,-2,-14,8,1,-7,1,-8,-5,-3,-4,-13,-3,-5,-5,-13,-12,-7],
             [-2,-9,-2,2,-14,1,8,-4,-5,-5,-9,-4,-7,-14,-5,-4,-6,-17,-8,-6],
             [-2,-9,-3,-3,-9,-7,-4,6,-9,-11,-10,-7,-8,-9,-6,-2,-6,-15,-14,-5],
             [-7,-2,0,-4,-7,1,-5,-9,9,-9,-6,-6,-10,-6,-4,-6,-7,-7,-3,-6],
             [-5,-5,-5,-7,-6,-8,-5,-11,-9,8,-1,-6,-1,-2,-8,-7,-2,-14,-6,2],
             [-6,-8,-7,-12,-15,-5,-9,-10,-6,-1,7,-8,1,-3,-7,-8,-7,-6,-7,-2],
             [-7,0,-1,-4,-14,-3,-4,-7,-6,-6,-8,7,-2,-14,-6,-4,-3,-12,-9,-9],
             [-5,-4,-9,-11,-13,-4,-7,-8,-10,-1,1,-2,11,-4,-8,-5,-4,-13,-11,-1],
             [-8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8],
             [-2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-7,-6,-8,-10,8,-2,-4,-14,-13,-6],
             [0,-3,0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6],
             [-1,-6,-2,-5,-8,-5,-6,-6,-7,-2,-7,-3,-4,-9,-4,0,7,-13,-6,-3],
             [-13,-2,-8,-15,-15,-13,-17,-15,-7,-14,-6,-12,-13,-4,-14,-5,-13,13,-5,-15],
             [-8,-10,-4,-11,-4,-12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7],
             [-2,-8,-8,-8,-6,-7,-6,-5,-6,2,-2,-9,-1,-8,-6,-6,-3,-15,-7,7])
    
    aa_map = {
        'A' : '0',
        'R' : '1',
        'N' : '2',
        'D' : '3',
        'C' : '4',
        'Q' : '5',
        'E' : '6',
        'G' : '7',
        'H' : '8',
        'I' : '9',
        'L' : '10',
        'K' : '11',
        'M' : '12',
        'F' : '13',
        'P' : '14',
        'S' : '15',
        'T' : '16',
        'W' : '17',
        'Y' : '18',
        'V' : '19'
        }
    score = 0
    for x in xrange(0, len(cdr_loop_seq)):
        score += PAM30[int(aa_map[cdr_loop_seq[x:x+1]])][int(aa_map[template_loop_seq[x:x+1]])]
    return score


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
myPath= os.getcwd()

profit_infile = 'profit.in'
create_profit_infile(profit_infile)


'''
cdrpdbtypes = ['A_cdr1', 'A_cdr2', 'A_cdr3', 'A_hv4', 'B_cdr1', 'B_cdr2', 'B_cdr3', 'B_hv4']
cdrlentypes = range(1, 26)
for f in cdrpdbtypes:
    for t in cdrlentypes:
        ftag = "*_*_"+str(f)+"_"+str(t)+".pdb"
        if len(glob.glob1(myPath,ftag)) > 0:
            print ftag, len(glob.glob1(myPath,ftag))
            ofile = str(f)+"_"+str(t)+".rms"
            gen_pairwise_rmsd_matrix(myPath,ftag,ofile)
            

'''
cdrlentypes = range(11, 13)
for f in cdrlentypes:   
    tcr_types = ('*_B_cdr3_'+str(f)+'.pdb')
    ab_types = ('*_H3_'+str(f)+'.pdb')
    tcr_files_grabbed = []
    tcr_files_grabbed.extend(glob.glob(tcr_types))
    ab_files_grabbed = []
    ab_files_grabbed.extend(glob.glob(ab_types))
    if len(tcr_files_grabbed) > 0:
        if len(ab_files_grabbed) > 0:
            for i in tcr_files_grabbed:
                inpose1 = rosetta.Pose()
                rosetta.core.import_pose.pose_from_pdb( inpose1 , i)
                inpose1_seq = inpose1.sequence()
                for j in ab_files_grabbed:
                    rmsd = get_rms(profit_infile, i, j)
                    if (float(rmsd) < 1.0):
                        inpose2 = rosetta.Pose()
                        rosetta.core.import_pose.pose_from_pdb( inpose2 , j)
                        score = calc_alignment_score(inpose1.sequence(), inpose2.sequence())
                        print str(rmsd).rstrip(), str(score).rstrip(), i, j
