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
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocols -ignore_zero_occupancy false")
import rosetta.protocols.grafting as graft


import gzip

script, fasta_file, TAG = argv

def extract_tcr_domain(inpose, TAG, record):
    seq = inpose.sequence()
    process = Popen(["ANARCI", "-i", seq, "-o", "tmp.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    anarcilist = []
    with open("tmp.out") as f:
        for line in f:
            if line.startswith("@"):
                #print line
                values = line.split("|")
                if values[3] == TAG:
                    seqstart = values[6]
                    seqend = values[7]
            if line.startswith(TAG):
                values = [line[0],line[2:7],line[8],line[10]]
                if line[10] == "-":continue
                anarcilist.append(values)

    domainpose = rosetta.Pose()
    domainpose = graft.return_region( inpose, int(seqstart)+1, int(seqend)+1)
    if domainpose.total_residue() == len(anarcilist):
        domainpose.dump_pdb("temp1.pdb")    
        OUT = open("temp2.pdb", 'w+')
        with open("temp1.pdb") as IN:
            for line in IN.readlines():
                if line[0:4] == 'ATOM':
                    slnum =  int(line[22:26].strip())
                    slnum -= 1
                    line  = line[:22] + anarcilist[slnum][1].strip().rjust(4) + anarcilist[slnum][2] + line[27:]
                    OUT.write(line)
    
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
 
def extract_chain_and_seq_from_structure(structure, chainid, out_file):

    #to remove HETATM records
    class nonHetSelect(Select):
        def accept_residue(self,residue):
            if residue.id[0] == ' ':
                return 1
            else:
                return 0

    pdbfile = parser.get_structure("PDB", structure)
    mychain = pdbfile[0][chainid]
    io.set_structure(mychain)
    io.save(out_file+'.pdb', nonHetSelect())

def find_cdr(seq, TAG):
    seq = str(seq.rstrip())
    process = Popen(["ANARCI", "-i", seq, "-o", "res.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)

    ahonumlist = []
    count = 0
    #alpha
    cdr1_begin_pos = None
    cdr2_begin_pos = None
    cdr3_begin_pos = None
    hv4_begin_pos = None
    cdr1_end_pos = None
    cdr2_end_pos = None
    cdr3_end_pos = None
    hv4_end_pos = None

    with open("res.out") as f:
        for line in f:
            if line.startswith("@"):
                values = line.split("|")
                if values[3]  == TAG:
                    seqstart = values[6]
                    seqend = values[7]
            #if not ( line.startswith("#") or line.startswith("@") or line.startswith("//") )
            if line.startswith(TAG):
                values = [line[0],line[2:7],line[8],line[10]]
                if values[3] == "-":continue
                if int(values[1].strip()) == cdr_aho_num[0]: cdr1_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[1]: cdr1_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[2]: cdr2_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[3]: cdr2_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[4]: cdr3_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[5]: cdr3_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[6]: hv4_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[7]: hv4_end_pos = count
               
                count += 1                
                ahonumlist.append(values)

    cdr1 = ''
    cdr2 = ''
    cdr3 = ''
    hv4 = ''


    for item in ahonumlist[cdr1_begin_pos:cdr1_end_pos+1]:
        cdr1 += item[3].rstrip()
    for item in ahonumlist[cdr2_begin_pos:cdr2_end_pos+1]:
        cdr2 += item[3].rstrip()
    for item in ahonumlist[cdr3_begin_pos:cdr3_end_pos+1]:
        cdr3 += item[3].rstrip()
    for item in ahonumlist[hv4_begin_pos:hv4_end_pos+1]:
        hv4 += item[3].rstrip()
    
    return cdr1,cdr2,cdr3,hv4,seqstart,seqend,cdr1_begin_pos,cdr1_end_pos,cdr2_begin_pos,cdr2_end_pos,cdr3_begin_pos,cdr3_end_pos,hv4_begin_pos,hv4_end_pos

#####
#MAIN
#####

if TAG == "A":
    cdr_aho_num = [24,42,56,71,107,138,81,89]
elif TAG == "B":
    cdr_aho_num = [24,42,56,70,107,138,81,89]
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    #seq based
    print "PDB", record.id
    cdrseq = find_cdr(record.seq,TAG)

    #str based
    #Extract chain from PDB
    gzpdbfile_path = '/TCRmodeller/PDB_RELEASE/pdb_structures' + '/%s/pdb%s.ent.gz' %(record.id[1:3], record.id[:4])
    gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
    extract_chain_and_seq_from_structure(gzpdbfile, record.id[5:6], record.id[:4]+"_"+record.id[5:6])
    
    inpose = rosetta.Pose()
    infile = record.id[:4]+"_"+record.id[5:6]+".pdb"
    rosetta.core.import_pose.pose_from_pdb( inpose , infile )
    extract_tcr_domain(inpose,TAG,record)
    ahopose = rosetta.Pose()
    rosetta.core.import_pose.pose_from_pdb( ahopose , "temp2.pdb" )

    chainid = "A"#rosetta output pose has default chain id as 'A'
    if TAG == "A": chaintype = 'Alpha'
    elif TAG == "B": chaintype = 'Beta'
    else: chaintype = None
    
    res = extract_region(ahopose, cdr_aho_num[0], cdr_aho_num[1], chainid)    
    if res == "ERROR":
        print "cdr1",cdrseq[0],record.id[:4],record.id[5:6],len(cdrseq[0]),chaintype,"Mismatch"
    elif res.total_residue() != len(cdrseq[0]):
        print "cdr1",cdrseq[0],record.id[:4],record.id[5:6],len(cdrseq[0]),chaintype,"Mismatch"
    else:
        #print cdrseq[0], res.sequence()
        assert (cdrseq[0] == res.sequence())
        print "cdr1",cdrseq[0],record.id[:4],record.id[5:6],len(cdrseq[0]),chaintype,"OK"        
        outfile = record.id[:4]+"_"+record.id[5:6]+"_"+TAG+"_cdr1"+"_"+str(len(cdrseq[0]))+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(ahopose, cdr_aho_num[2], cdr_aho_num[3], chainid)    
    if res == "ERROR":
        print "cdr2",cdrseq[1],record.id[:4],record.id[5:6],len(cdrseq[1]),chaintype,"Mismatch"
    elif res.total_residue() != len(cdrseq[1]):
        print "cdr2",cdrseq[1],record.id[:4],record.id[5:6],len(cdrseq[1]),chaintype,"Mismatch"
    else:
        #print cdrseq[1], res.sequence()
        assert (cdrseq[1] == res.sequence())
        print "cdr2",cdrseq[1],record.id[:4],record.id[5:6],len(cdrseq[1]),chaintype,"OK"        
        outfile = record.id[:4]+"_"+record.id[5:6]+"_"+TAG+"_cdr2"+"_"+str(len(cdrseq[1]))+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(ahopose, cdr_aho_num[4], cdr_aho_num[5], chainid)    
    if res == "ERROR":
        print "cdr3",cdrseq[2],record.id[:4],record.id[5:6],len(cdrseq[2]),chaintype,"Mismatch"
    elif res.total_residue() != len(cdrseq[2]):
        print "cdr3",cdrseq[2],record.id[:4],record.id[5:6],len(cdrseq[2]),chaintype,"Mismatch"
    else:
        #print cdrseq[2], res.sequence()
        assert (cdrseq[2] == res.sequence())
        print "cdr3",cdrseq[2],record.id[:4],record.id[5:6],len(cdrseq[2]),chaintype,"OK"        
        outfile = record.id[:4]+"_"+record.id[5:6]+"_"+TAG+"_cdr3"+"_"+str(len(cdrseq[2]))+".pdb"
        res.dump_pdb(outfile)

    res = extract_region(ahopose, cdr_aho_num[6], cdr_aho_num[7], chainid)    
    if res == "ERROR":
        print "hv4",cdrseq[3],record.id[:4],record.id[5:6],len(cdrseq[3]),chaintype,"Mismatch"
    elif res.total_residue() != len(cdrseq[3]):
        print "hv4",cdrseq[3],record.id[:4],record.id[5:6],len(cdrseq[3]),chaintype,"Mismatch"
    else:
        #print cdrseq[3], res.sequence()
        assert (cdrseq[3] == res.sequence())
        print "hv4",cdrseq[3],record.id[:4],record.id[5:6],len(cdrseq[3]),chaintype,"OK"        
        outfile = record.id[:4]+"_"+record.id[5:6]+"_"+TAG+"_hv4"+"_"+str(len(cdrseq[3]))+".pdb"
        res.dump_pdb(outfile)


