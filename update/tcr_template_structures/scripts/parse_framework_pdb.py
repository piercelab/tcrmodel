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

#####
#MAIN
#####

cdra_aho_num = [24,42,57,71,107,138,81,89]
cdrb_aho_num = [24,42,57,70,107,138,81,89]

with open(inp_file) as f:
    for infile in f:
#        print infile
        infile = infile.strip()
        infile_base = os.path.splitext(infile)[0]
        #str based
        inpose = rosetta.Pose()
        rosetta.core.import_pose.pose_from_pdb( inpose , infile )
        seq = inpose.sequence()
        #print seq
        process = Popen(["ANARCI", "-i", seq, "-o", "tmp1.out", "-r", "tr", "-s" , "a"  ], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        #print(stdout)
        #print(stderr)
        regexa = "\A#\|.*\|"+"A"+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
        anarcilista = []
        counta = 0
        start_posa = ''
        end_posa = ''


        regexb = "\A#\|.*\|"+"B"+"\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
        anarcilistb = []
        countb = 0
        start_posb = ''
        end_posb = ''

    
        with open("tmp1.out") as f:
            for line in f:
                resa = re.search(regexa, line)
                if resa:
                    seqstarta = resa.groups()[0]
                    seqenda = resa.groups()[1]

                resb = re.search(regexb, line)
                if resb:
                    seqstartb = resb.groups()[0]
                    seqendb = resb.groups()[1]

                if line.startswith("A"):
                    if start_posa == '': start_posa = line[2:7].strip()
                    end_posa = line[2:7].strip()
                    if line[10] == "-":continue
                    valuesa = [line[0],line[2:7],line[8],line[10]]                

                    if line[2:7].strip() == "56": posa_56 = counta
                    if line[2:7].strip() == "90": posa_90 = counta
                    if int(line[2:7].strip()) == cdra_aho_num[0]: cdr1a_begin_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[1]: cdr1a_end_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[2]: cdr2a_begin_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[3]: cdr2a_end_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[4]: cdr3a_begin_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[5]: cdr3a_end_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[6]: hv4a_begin_pos = counta
                    if int(line[2:7].strip()) == cdra_aho_num[7]: hv4a_end_pos = counta
                    anarcilista.append(valuesa)
                    counta += 1

                if line.startswith("B"):
                    if start_posb == '': start_posb = line[2:7].strip()
                    end_posb = line[2:7].strip()
                    if line[10] == "-":continue
                    valuesb = [line[0],line[2:7],line[8],line[10]]                

                    if line[2:7].strip() == "56": posb_56 = countb
                    if line[2:7].strip() == "90": posb_90 = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[0]: cdr1b_begin_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[1]: cdr1b_end_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[2]: cdr2b_begin_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[3]: cdr2b_end_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[4]: cdr3b_begin_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[5]: cdr3b_end_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[6]: hv4b_begin_pos = countb
                    if int(line[2:7].strip()) == cdrb_aho_num[7]: hv4b_end_pos = countb
                    anarcilistb.append(valuesb)
                    countb += 1

        num_begin = 19
        num_end = 9
        start_cap = 0
        end_cap = 0
        if (len(str_vdomain[:cdr1_begin_pos]) > num_begin):
            start_cap = len(str_vdomain[:cdr1_begin_pos]) - num_begin
            anarcilist = anarcilist[start_cap:]
        elif (len(str_vdomain[:cdr1_begin_pos]) == num_begin):
            pass
        else:
            continue
        
        if (len(str_vdomain[cdr3_end_pos+1:]) > num_end):
            end_cap = len(str_vdomain[cdr3_end_pos+1:]) - num_end
            anarcilist = anarcilist[:-end_cap]
        elif(len(str_vdomain[cdr3_end_pos+1:]) == num_end):
            pass#anarcilist = anarcilist                                          
        else:
            continue
        
        domainposea = rosetta.Pose()
        domainposea = graft.return_region( inpose, int(seqstarta)+1+int(posa_3), int(seqenda)+1)
        domainposea.dump_pdb(infile_base+"_Adomain.pdb")    
        str_Adomain = domainposea.sequence()
            
        OUT = open(infile_base+"_Aaho.pdb", 'w+')
        with open(infile_base+"_Adomain.pdb") as IN:
            for pdbline in IN.readlines():
                if pdbline[0:4] == 'ATOM':
                    slnum =  int(pdbline[22:26].strip())
                    slnum -= 1
                    pdbline  = pdbline[:22] + anarcilista[slnum][1].strip().rjust(4) + anarcilista[slnum][2] + pdbline[27:]
                    OUT.write(pdbline)

        domainposeb = rosetta.Pose()
        domainposeb = graft.return_region( inpose, int(seqstartb)+1+int(posb_3), int(seqendb)+1)
        domainposeb.dump_pdb(infile_base+"_Bdomain.pdb")    
        str_Bdomain = domainposeb.sequence()
            #print "str_Bdomain" , str_Bdomain
        if not posb_3 == 0:
            anarcilistb = anarcilistb[posb_3:]
            
            OUT = open(infile_base+"_Baho.pdb", 'w+')
            with open(infile_base+"_Bdomain.pdb") as IN:
                for pdbline in IN.readlines():
                    if pdbline[0:4] == 'ATOM':
                        slnum =  int(pdbline[22:26].strip())
                        slnum -= 1
                        pdbline  = pdbline[:22] + anarcilistb[slnum][1].strip().rjust(4) + anarcilistb[slnum][2] + pdbline[27:]
                        OUT.write(pdbline)


            print str_Adomain
            print str_Adomain[:cdr1a_begin_pos], str_Adomain[cdr1a_begin_pos:cdr1a_end_pos+1], str_Adomain[cdr1a_end_pos+1:posa_56], str_Adomain[posa_56:posa_90+1], str_Adomain[posa_90+1:cdr3a_begin_pos], str_Adomain[cdr3a_begin_pos:cdr3a_end_pos+1], str_Adomain[cdr3a_end_pos+1:]

            print str_Bdomain
            print str_Bdomain[:cdr1b_begin_pos], str_Bdomain[cdr1b_begin_pos:cdr1b_end_pos+1], str_Bdomain[cdr1b_end_pos+1:posb_56], str_Bdomain[posb_56:posb_90+1], str_Bdomain[posb_90+1:cdr3b_begin_pos], str_Bdomain[cdr3b_begin_pos:cdr3b_end_pos+1], str_Bdomain[cdr3b_end_pos+1:]

#            frameworkseqa = str_Adomain[:cdr1a_begin_pos] + str_Adomain[cdr1a_end_pos+1:posa_56] + str_Adomain[posa_90+1:cdr3a_begin_pos] + str_Adomain[cdr3a_end_pos+1:]
#            frameworkseqb = str_Bdomain[:cdr1b_begin_pos] + str_Bdomain[cdr1b_end_pos+1:posb_56] + str_Bdomain[posb_90+1:cdr3b_begin_pos] + str_Bdomain[cdr3b_end_pos+1:]
            frameworkseqa = str_Adomain[cdr1a_end_pos+1:posa_56] + str_Adomain[posa_90+1:cdr3a_begin_pos]
            frameworkseqb = str_Bdomain[cdr1b_end_pos+1:posb_56] + str_Bdomain[posb_90+1:cdr3b_begin_pos]


            with open('fw_orientation.fasta', 'a') as fr:
                fr.write(infile_base+" "+frameworkseqa+" "+frameworkseqb+"\n")
                print infile_base+" "+frameworkseqa+" "+frameworkseqb

                

