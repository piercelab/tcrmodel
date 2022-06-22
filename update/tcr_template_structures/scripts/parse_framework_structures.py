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
#rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -ignore_zero_occupancy false")
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocols -ignore_zero_occupancy false")
import rosetta.protocols.grafting as graft

import gzip

script, inp_file = argv

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

def parse_tcr_seq(seq,recordid,TAG):
    seq = str(seq.rstrip())
    process = Popen(["ANARCI", "-i", seq, "-o", "res.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    
    start_pos = ''
    end_pos = ''
    count = 0
    regex0 = "\A#\|.*\|"+TAG+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    with open("res.out") as f:
        for line in f:
            res = re.search(regex0, line)
            if res:
                tcrdomain = seq[int(res.groups()[0]):int(res.groups()[1])+1]
                with open(TAG+'_tcr_vdomain.fasta', 'a') as f:
                    f.write(">"+recordid+"\n")
                    f.write(tcrdomain+"\n")
                    
            if line.startswith(TAG):
                if start_pos == '': start_pos = line[2:7].strip()
                end_pos = line[2:7].strip()
                
                if line[10] == "-":continue
                
                if line[2:7].strip() == "55": pos_55 = count
                if line[2:7].strip() == "91": pos_91 = count
                
                if int(line[2:7].strip()) == cdr_aho_num[0]: cdr1_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[1]: cdr1_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[2]: cdr2_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[3]: cdr2_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[4]: cdr3_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[5]: cdr3_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[6]: hv4_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[7]: hv4_end_pos = count
                count += 1

#    print "domain", tcrdomain    
#    print "CDRs", tcrdomain[cdr1_begin_pos:cdr1_end_pos+1],  tcrdomain[cdr2_begin_pos:cdr2_end_pos+1],  tcrdomain[cdr3_begin_pos:cdr3_end_pos+1] 
    frameworkseq = tcrdomain[:cdr1_begin_pos]+tcrdomain[cdr1_end_pos+1:cdr2_begin_pos]+tcrdomain[cdr2_end_pos+1:cdr3_begin_pos]+tcrdomain[cdr3_end_pos+1:]
#    print "start_pos", start_pos, recordid
    frameworkseq = tcrdomain[:cdr1_begin_pos]+tcrdomain[cdr1_end_pos+1:pos_55+1]+tcrdomain[pos_91:cdr3_begin_pos]+tcrdomain[cdr3_end_pos+1:]


    cdr3hv4 = tcrdomain[pos_55+1:pos_91]
#    print "cdr3hv4", cdr3hv4
    with open(TAG+'_tcr_cdr3hv4.fasta', 'a') as fr:
        fr.write(">"+recordid+"\n")
        fr.write(cdr3hv4+"\n")
    
    return tcrdomain, frameworkseq, start_pos, end_pos

def extract_framework(seq, TAG, outtag):      

    if TAG == "A":
        cdr_aho_num = [24,42,57,71,107,138,81,89]
    elif TAG == "B":
        cdr_aho_num = [24,42,57,70,107,138,81,89]

    process = Popen(["ANARCI", "-i", seq, "-o", "tmp1.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    anarcilist = []
    count = 0
    start_pos = ''
    end_pos = ''
    regex0 = "\A#\|.*\|"+TAG+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"

    with open("tmp1.out") as f:
        for line in f:
            res = re.search(regex0, line)
            if res:
                seqstart = res.groups()[0]
                seqend = res.groups()[1]

            if line.startswith(TAG):
                if start_pos == '': start_posa = line[2:7].strip()
                end_pos = line[2:7].strip()
                if line[10] == "-":continue
                values = [line[0],line[2:7],line[8],line[10]]
            
                if line[2:7].strip() == "56": pos_56 = count
                if line[2:7].strip() == "90": pos_90 = count
                if int(line[2:7].strip()) == cdr_aho_num[0]: cdr1_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[1]: cdr1_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[2]: cdr2_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[3]: cdr2_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[4]: cdr3_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[5]: cdr3_end_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[6]: hv4_begin_pos = count
                if int(line[2:7].strip()) == cdr_aho_num[7]: hv4_end_pos = count
                anarcilist.append(values)
                count += 1

    domainpose = rosetta.Pose()
    domainpose = graft.return_region( inpose, int(seqstart)+1, int(seqend)+1)
    domainpose.dump_pdb(outtag+"_vdomain.pdb")    
    str_vdomain = domainpose.sequence()

    num_begin = 20
    num_end = 10
    start_cap = 0
    end_cap = 0
    if (len(str_vdomain[:cdr1_begin_pos]) > num_begin):
        start_cap = len(str_vdomain[:cdr1_begin_pos]) - num_begin
        anarcilist = anarcilist[start_cap:]
    elif (len(str_vdomain[:cdr1_begin_pos]) == num_begin):
        pass
    else:
        return None
    
    if (len(str_vdomain[cdr3_end_pos+1:]) > num_end):
        end_cap = len(str_vdomain[cdr3_end_pos+1:]) - num_end
        anarcilist = anarcilist[:-end_cap]
    elif(len(str_vdomain[cdr3_end_pos+1:]) == num_end):
        pass#anarcilist = anarcilist
    else:
        return None

    newdomainpose = rosetta.Pose()
    newdomainpose = graft.return_region( inpose, int(seqstart)+1+start_cap, int(seqend)+1-end_cap)
    newdomainpose.dump_pdb(outtag+"_"+TAG+".pdb")

    frameworkseq = str_vdomain[start_cap:cdr1_begin_pos] + str_vdomain[cdr1_end_pos+1:pos_56] + str_vdomain[pos_90+1:cdr3_begin_pos] + str_vdomain[cdr3_end_pos+1:cdr3_end_pos+1+num_end]
    
    OUT = open(outtag+"_"+TAG+"_aho.pdb", 'w+')
    with open(outtag+"_"+TAG+".pdb") as IN:
        for pdbline in IN.readlines():
            if pdbline[0:4] == 'ATOM':
                slnum =  int(pdbline[22:26].strip())
                slnum -= 1
                pdbline  = pdbline[:22] + anarcilist[slnum][1].strip().rjust(4) + anarcilist[slnum][2] + pdbline[27:]
                OUT.write(pdbline)

    return frameworkseq


#####
#MAIN
#####

with open(inp_file) as f:
    for infile in f:
        infile = infile.strip()
        infile_base = os.path.splitext(infile)[0]
        
        inpose = rosetta.Pose()
        rosetta.core.import_pose.pose_from_pdb( inpose , infile )
        seq = inpose.sequence()

        fwseqa = extract_framework(seq, "A", infile_base)
        fwseqb = extract_framework(seq, "B", infile_base)
        if fwseqa and fwseqb is not None:
            print infile_base+" "+fwseqa+" "+fwseqb
            with open('tcr_framework.seq', 'a') as fr:
                fr.write(infile_base+" "+fwseqa+" "+fwseqb+"\n")
        else:
            print infile_base, "None"
