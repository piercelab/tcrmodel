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

script, fasta_file, TAG = argv

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
    process = Popen(["ANARCI", "-i", seq, "-o", "anr_seq.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    
    start_pos = ''
    end_pos = ''
    count = 0
    regex0 = "\A#\|.*\|"+TAG+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    with open("anr_seq.out") as f:
        for line in f:
            res = re.search(regex0, line)
            if res:
                tcrdomain = seq[int(res.groups()[0]):int(res.groups()[1])+1]
                with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as f:
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

    frameworkseq = tcrdomain[:cdr1_begin_pos]+tcrdomain[cdr1_end_pos+1:pos_55+1]+tcrdomain[pos_91:cdr3_begin_pos]+tcrdomain[cdr3_end_pos+1:]

    cdr1 = tcrdomain[cdr1_begin_pos:cdr1_end_pos+1]
    cdr2hv4 = tcrdomain[pos_55+1:pos_91]
    cdr3 = tcrdomain[cdr3_begin_pos:cdr3_end_pos+1]

    return tcrdomain, frameworkseq, cdr1, cdr2hv4, cdr3

#####
#MAIN
#####

if TAG == "A":
    cdr_aho_num = [24,42,57,71,107,138,81,89]
elif TAG == "B":
    cdr_aho_num = [24,42,57,70,107,138,81,89]
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    outtag = record.id[:4]+"_"+record.id[5:6]
    print outtag
    
    #seq based
    seq_res = parse_tcr_seq(record.seq,record.id,TAG)    
    seq_vdomain = seq_res[0]    
    seq_framework = seq_res[1]    
    seq_cdr1 = seq_res[2]    
    seq_cdr2hv4 = seq_res[3]    
    seq_cdr3 = seq_res[4]    

    #str based
    #Extract chain from PDB
    gzpdbfile_path = '/TCRmodeller/PDB_RELEASE/pdb_structures' + '/%s/pdb%s.ent.gz' %(record.id[1:3], record.id[:4])
    gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
    extract_chain_and_seq_from_structure(gzpdbfile, record.id[5:6], outtag)

    inpose = rosetta.Pose()
    infile = outtag+".pdb"
    rosetta.core.import_pose.pose_from_pdb( inpose , infile )
    seq = inpose.sequence()
    process = Popen(["ANARCI", "-i", seq, "-o", "anr_str.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #print(stdout)
    #print(stderr)
    anarcilist = []
    count = 0
    start_pos = ''
    end_pos = ''
    regex0 = "\A#\|.*\|"+TAG+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    seqstart = None
    seqend = None
    with open("anr_str.out") as f:
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

    if seqstart is None and seqend is None:
        print outtag, "ANARCI failed for the sequence generated from PDB file"            
        continue

    #output vdomain pdb
    domainpose = rosetta.Pose()
    domainpose = graft.return_region( inpose, int(seqstart)+1, int(seqend)+1)
    domainpose.dump_pdb(outtag+"_vdomain.pdb")    
    str_vdomain = domainpose.sequence()
    str_cdr1 = str_vdomain[cdr1_begin_pos:cdr1_end_pos+1]
    str_cdr3 = str_vdomain[cdr3_begin_pos:cdr3_end_pos+1]
    str_cdr2hv4 = str_vdomain[pos_56:pos_90+1]
    str_framework = str_vdomain[:cdr1_begin_pos]+str_vdomain[cdr1_end_pos+1:pos_56] + str_vdomain[pos_90+1:cdr3_begin_pos]+str_vdomain[cdr3_end_pos+1:]
    print seq_framework, seq_cdr1, seq_cdr2hv4, seq_cdr3
    print str_framework, str_cdr1, str_cdr2hv4, str_cdr3

    #output cdr1,cdr2hv4 and cdr3 fasta
    cdr1_check = None
    cdr2hv4_check = None
    cdr3_check = None
    if(seq_cdr1 == str_cdr1):
        cdr1_check = True        
        with open(TAG+'_TCR_CDR1_'+str(len(str_cdr1))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"\n")
            fr.write(str_cdr1+"\n")
    if(seq_cdr2hv4 == str_cdr2hv4):
        cdr2hv4_check = True
        with open(TAG+'_TCR_CDR2HV4_'+str(len(str_cdr2hv4))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"\n")
            fr.write(str_cdr2hv4+"\n")
    if(seq_cdr3 == str_cdr3):
        cdr3_check = True
        with open(TAG+'_TCR_CDR3_'+str(len(str_cdr3))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"\n")
            fr.write(str_cdr3+"\n")

    #cap the begin and end of framework sequence
    num_begin = 20
    num_end = 10

    #start cap
    start_cap = 0
    startcap_check = True
    if (len(str_vdomain[:cdr1_begin_pos]) > num_begin):
        start_cap = len(str_vdomain[:cdr1_begin_pos]) - num_begin
        capped_anarcilist = anarcilist[start_cap:]
    elif (len(str_vdomain[:cdr1_begin_pos]) == num_begin):
        pass
    else:
        print "start cap failed:", num_begin, str_vdomain[:cdr1_begin_pos], len(str_vdomain[:cdr1_begin_pos])
        startcap_check = False

    #end_Cap
    end_cap = 0
    endcap_check = True;
    if (len(str_vdomain[cdr3_end_pos+1:]) > num_end):
        end_cap = len(str_vdomain[cdr3_end_pos+1:]) - num_end
        capped_anarcilist = anarcilist[:-end_cap]
    elif(len(str_vdomain[cdr3_end_pos+1:]) == num_end):
        pass
    else:
        print "end_cap_failed:", num_end, str_vdomain[cdr3_end_pos+1:], len(str_vdomain[cdr3_end_pos+1:])
        endcap_check = False

    #output capped framework sequence
    if (startcap_check and endcap_check):
        if((seq_vdomain == str_vdomain) or (re.search(str_vdomain, seq_vdomain) is not None)):
            framework_check = True
            frameworkseq = str_vdomain[start_cap:cdr1_begin_pos] + str_vdomain[cdr1_end_pos+1:pos_56] + str_vdomain[pos_90+1:cdr3_begin_pos] + str_vdomain[cdr3_end_pos+1:cdr3_end_pos+1+num_end]
            print "framework_check", framework_check, frameworkseq
            #with open(TAG+'_TCR_FW_'+str(len(frameworkseq))+'.fasta', 'a') as fr:
            with open(TAG+'_TCR_FW.fasta', 'a') as fr:
                fr.write(">"+record.id+"\n")
                fr.write(frameworkseq+"\n")

        #output capped aho pdb
            newdomainpose = rosetta.Pose()
            newdomainpose = graft.return_region( inpose, int(seqstart)+1+start_cap, int(seqend)+1-end_cap)
            newdomainpose.dump_pdb(outtag+"_capped_"+TAG+".pdb")
            OUT = open(outtag+"_capped_"+TAG+"_aho.pdb", 'w+')
            with open(outtag+"_capped_"+TAG+".pdb") as IN:
                for pdbline in IN.readlines():
                    if pdbline[0:4] == 'ATOM':
                        slnum =  int(pdbline[22:26].strip())
                        slnum -= 1
                        pdbline  = pdbline[:22] + capped_anarcilist[slnum][1].strip().rjust(4) + capped_anarcilist[slnum][2] + pdbline[27:]
                        OUT.write(pdbline)
        else:     
            print "seq_vdomain and str_vdomain are not same"
            print "seq_vdomain", seq_vdomain
            print "str_vdomain", str_vdomain
            framework_check = False

    #output vdomain aho pdb
    print "check : ", framework_check, cdr1_check, cdr2hv4_check, cdr3_check
    if (framework_check or cdr1_check or cdr2hv4_check or cdr3_check):
        OUT = open(outtag+"_vdomain_aho.pdb", 'w+')
        with open(outtag+"_vdomain.pdb") as IN:
            for pdbline in IN.readlines():
                if pdbline[0:4] == 'ATOM':
                    slnum =  int(pdbline[22:26].strip())
                    slnum -= 1
                    pdbline  = pdbline[:22] + anarcilist[slnum][1].strip().rjust(4) + anarcilist[slnum][2] + pdbline[27:]
                    OUT.write(pdbline)

