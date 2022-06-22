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

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocols -ignore_zero_occupancy false -ignore_unrecognized_res")
import rosetta.protocols.grafting as graft

import glob

import gzip

#script, pdblist, TAG = argv

    
def extract_region(ahopose, pdb_start_pos, pdb_end_pos, chainid):
    pose_start_pos = ahopose.pdb_info().pdb2pose(chainid, pdb_start_pos)
    pose_end_pos = ahopose.pdb_info().pdb2pose(chainid, pdb_end_pos)
    #print pdb_start_pos, pdb_end_pos, pose_start_pos, pose_end_pos
    if (pose_start_pos == 0) or (pose_end_pos == 0):
        return "ERROR"
    else:
        cdrpose = rosetta.Pose()
        cdrpose = graft.return_region( ahopose, pose_start_pos, pose_end_pos)
        return cdrpose

#####
#MAIN
#####


L_cdr_aho_num = [24,42,57,71,107,138]
H_cdr_aho_num = [24,42,57,70,107,138]
    
    
myPath = "/home/rg/Fv/splitchains/"
listfiles = glob.glob1(myPath,"*_L_aho.pdb")

for infile in listfiles:
    
    print infile    
    inpose = rosetta.Pose()
    rosetta.core.import_pose.pose_from_pdb( inpose , myPath+infile )
    
    bname = os.path.splitext(os.path.basename(infile))[0]
    
    res = extract_region(inpose, L_cdr_aho_num[4], L_cdr_aho_num[5], "L")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_L3"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)
    


'''
    res = extract_region(inpose, L_cdr_aho_num[0], L_cdr_aho_num[1], "L")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_L1"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)
        
    res = extract_region(inpose, L_cdr_aho_num[2], L_cdr_aho_num[3], "L")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_L2"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(inpose, H_cdr_aho_num[0], H_cdr_aho_num[1], "H")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_H1"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(inpose, H_cdr_aho_num[2], H_cdr_aho_num[3], "H")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_H2"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(inpose, H_cdr_aho_num[4], H_cdr_aho_num[5], "H")    
    if not res == "ERROR":
        cdrlen = res.total_residue()
        outfile = bname+"_H3"+"_"+str(cdrlen)+".pdb"
        res.dump_pdb(outfile)
'''
