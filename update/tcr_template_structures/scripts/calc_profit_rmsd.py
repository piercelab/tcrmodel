#!/usr/bin/env python                                                                                                                       
import sys, os
from sys import argv

import subprocess
from subprocess import Popen, PIPE

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocols -ignore_zero_occupancy false -renumber_pdb")

#from TCR_functions import *

def create_profit_infile(infile):
    
    f = open(infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ATOMS N,CA,C\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE *:*\n")
    f.write("FIT\n")
    f.close()
    
def get_rms(infile, reference, mobile):
    rmslist = []
    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    cmd = subprocess.Popen([profit_program, '-f', infile, '-h', reference, mobile], stdout=subprocess.PIPE)
    for line in cmd.stdout:
        #print line
        rmsd = None
        if "RMS:" in line:
            rmsd = line.split(' ')[-1]
            rmslist.append(rmsd.rstrip())
    return rmslist
        

script, f1, f2 = argv

#renumber_pdbfile_to_aho(f1, None, None, True)
#renumber_pdbfile_to_aho(f2, None, None, True)
#ref = str(os.path.basename(f1)+".aho.pdb")
#mob = str(os.path.basename(f2)+".aho.pdb")


refpose = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( refpose , f1 )
ref_chaina = refpose.pdb_info().chain( 1 );
ref_bchain_start = refpose.conformation().chain_begin(2);
ref_chainb = refpose.pdb_info().chain( ref_bchain_start );
ref_cdr1as = refpose.pdb_info().pdb2pose(ref_chaina,24)
ref_cdr1ae = refpose.pdb_info().pdb2pose(ref_chaina,42)
ref_cdr2as = refpose.pdb_info().pdb2pose(ref_chaina,56)
ref_cdr2ae = refpose.pdb_info().pdb2pose(ref_chaina,71)
ref_cdr3as = refpose.pdb_info().pdb2pose(ref_chaina,107)
ref_cdr3ae = refpose.pdb_info().pdb2pose(ref_chaina,138)

ref_cdr1bs = refpose.pdb_info().pdb2pose(ref_chainb,24)
ref_cdr1be = refpose.pdb_info().pdb2pose(ref_chainb,42)
ref_cdr2bs = refpose.pdb_info().pdb2pose(ref_chainb,56)
ref_cdr2be = refpose.pdb_info().pdb2pose(ref_chainb,70)
ref_cdr3bs = refpose.pdb_info().pdb2pose(ref_chainb,107)
ref_cdr3be = refpose.pdb_info().pdb2pose(ref_chainb,138)

mobpose = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( mobpose , f2 )
mob_chaina = mobpose.pdb_info().chain( 1 );
mob_bchain_start = mobpose.conformation().chain_begin(2);
mob_chainb = mobpose.pdb_info().chain( mob_bchain_start );

mob_cdr1as = ref_cdr1as
mob_cdr1ae = ref_cdr1ae
mob_cdr2as = ref_cdr2as
mob_cdr2ae = ref_cdr2ae
mob_cdr3as = ref_cdr3as
mob_cdr3ae = ref_cdr3ae
mob_cdr1bs = ref_cdr1bs
mob_cdr1be = ref_cdr1be
mob_cdr2bs = ref_cdr2bs
mob_cdr2be = ref_cdr2be
mob_cdr3bs = ref_cdr3bs
mob_cdr3be = ref_cdr3be


ref_lena =  len(refpose.chain_sequence(1))
mob_lena =  len(mobpose.chain_sequence(1))
ref_tot =  refpose.total_residue()
mob_tot =  mobpose.total_residue()


ref_chaina = 'A'
ref_chainb = 'B'
mob_chaina = 'A'
mob_chainb = 'B'

infile ="profit.in"
f = open(infile,"w")
#f.write("ATOMS N,CA,C,O\n")
f.write("ATOMS N,CA,C\n")
f.write("FIT\n")
f.write("ZONE CLEAR\n")
f.write("ZONE "+ref_chaina+"1-"+ref_chaina+str(ref_cdr1as)+":"+mob_chaina+"1-"+mob_chaina+str(mob_cdr1as)+"\n")
f.write("ZONE "+ref_chaina+str(ref_cdr1ae)+"-"+ref_chaina+str(ref_cdr2as)+":"+mob_chaina+str(mob_cdr1ae)+"-"+mob_chaina+str(mob_cdr2as)+"\n")
f.write("ZONE "+ref_chaina+str(ref_cdr2ae)+"-"+ref_chaina+str(ref_cdr3as)+":"+mob_chaina+str(mob_cdr2ae)+"-"+mob_chaina+str(mob_cdr3as)+"\n")
f.write("ZONE "+ref_chaina+str(ref_cdr3ae)+"-"+ref_chaina+str(ref_lena)+":"+mob_chaina+str(mob_cdr3ae)+"-"+mob_chaina+str(mob_lena)+"\n")
f.write("ZONE "+ref_chainb+str((ref_lena)+1)+"-"+ref_chainb+str(ref_cdr1bs)+":"+mob_chainb+str((mob_lena)+1)+"-"+mob_chainb+str(mob_cdr1bs)+"\n")
f.write("ZONE "+ref_chaina+str(ref_cdr1be)+"-"+ref_chainb+str(ref_cdr2bs)+":"+mob_chainb+str(mob_cdr1be)+"-"+mob_chainb+str(mob_cdr2bs)+"\n")
f.write("ZONE "+ref_chainb+str(ref_cdr2be)+"-"+ref_chainb+str(ref_cdr3bs)+":"+mob_chainb+str(mob_cdr2be)+"-"+mob_chainb+str(mob_cdr3bs)+"\n")
f.write("ZONE "+ref_chainb+str(ref_cdr3be)+"-"+ref_chainb+str(ref_tot)+":"+mob_chainb+str(mob_cdr3be)+"-"+mob_chainb+str(mob_tot)+"\n")
f.write("FIT\n")
f.write("ZONE CLEAR\n")
f.write("DELRZONE ALL\n")
f.write("RZONE A"+str(ref_cdr1as)+"-"+ref_chaina+str(ref_cdr1ae)+":"+mob_chaina+str(mob_cdr1as)+"-"+mob_chaina+str(mob_cdr1ae)+"\n")
f.write("DELRZONE ALL\n")
f.write("RZONE A"+str(ref_cdr2as)+"-"+ref_chaina+str(ref_cdr2ae)+":"+mob_chaina+str(mob_cdr2as)+"-"+mob_chaina+str(mob_cdr2ae)+"\n")
f.write("DELRZONE ALL\n")
f.write("RZONE A"+str(ref_cdr3as)+"-"+ref_chaina+str(ref_cdr3ae)+":"+mob_chaina+str(mob_cdr3as)+"-"+mob_chaina+str(mob_cdr3ae)+"\n")
f.write("DELRZONE ALL\n")
f.write("RZONE B"+str(ref_cdr1bs)+"-"+ref_chainb+str(ref_cdr1be)+":"+mob_chainb+str(mob_cdr1bs)+"-"+mob_chainb+str(mob_cdr1be)+"\n")
f.write("DELRZONE ALL\n")
f.write("RZONE B"+str(ref_cdr2bs)+"-"+ref_chainb+str(ref_cdr2be)+":"+mob_chainb+str(mob_cdr2bs)+"-"+mob_chainb+str(mob_cdr2be)+"\n")
f.write("DELRZONE ALL\n")
f.write("RZONE B"+str(ref_cdr3bs)+"-"+ref_chainb+str(ref_cdr3be)+":"+mob_chainb+str(mob_cdr3bs)+"-"+mob_chainb+str(mob_cdr3be)+"\n")


f.close()
refpose.dump_pdb("ref.pdb")
mobpose.dump_pdb("mob.pdb")
rms = get_rms(infile, "ref.pdb", "mob.pdb")
print "RMSD", f1, f2, rms[0], rms[1], rms[2], rms[3], rms[4], rms[5], rms[6], rms[7]
