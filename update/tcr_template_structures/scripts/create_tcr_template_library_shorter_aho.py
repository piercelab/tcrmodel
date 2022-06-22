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

from TCR_functions import *
script, fasta_file, TAG = argv
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    pdbid = record.id[:4]
    chainid = record.id[5:6]
    outtag = pdbid+"_"+chainid+"_"+TAG
    pdbfile = outtag+".pdb"
    cap_pdbfile = outtag+"_cap.pdb"
    cap_ahofile = outtag+"_aho.pdb"
    nocap_ahofile = outtag+"_aho.pdb"
    print outtag
    #print record.seq
    #seq based
    fseq = record.seq
    if check_tcr_alpha_beta_domain(record.seq):
        fseq= get_tcr_domain_seq(record.seq,TAG)
    print fseq 
    regexres = assign_CDRs_using_REGEX_shorter_aho(fseq, TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seq_cdr1 = regexres.groups()[2]
        seq_cdr2hv4 = regexres.groups()[4]
        seq_cdr3 = regexres.groups()[6]   
        print "cap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
        with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_vdomain+"\n")
    else:
        regexres = assign_CDRs_using_REGEX_shorter_aho(fseq, TAG, True)
        if regexres:
            seq_vdomain = "None"
            seq_fw = "None"
            seq_cdr1 = regexres.groups()[2]
            seq_cdr2hv4 = regexres.groups()[4]
            seq_cdr3 = regexres.groups()[6]
            #print "nocap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
        else:
            print outtag, "FAIL: RegEx failed to identify TCR domain"
            continue

    get_pdb(record.id[:4],record.id[5:6],pdbfile,"/TCRmodeller/PDB_RELEASE/pdb_structures")
    strseq = pdb_to_fasta(pdbfile)
    pdbregex = re.search(seq_vdomain,strseq)
    if pdbregex:
        capbegin = pdbregex.span()[0]
        capend = pdbregex.span()[1]
        #print capbegin, capend
        cap_pdb(pdbfile,capbegin,capend,cap_pdbfile)
        renumber_pdbfile_to_aho(cap_pdbfile,TAG,cap_ahofile,True)
        with open(TAG+'_TCR_VDOMAIN_cap.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_vdomain+"\n")
        with open(TAG+'_TCR_FW.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_fw+"\n")
        with open(TAG+'_TCR_CDR1_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_cdr1+"\n")
        with open(TAG+'_TCR_CDR2HV4_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_cdr2hv4+"\n")
        with open(TAG+'_TCR_CDR3_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
            fr.write(">"+record.id+"_"+TAG+"\n")
            fr.write(seq_cdr3+"\n")
    else:
        tmp_ahofile = "temp.aho.pdb"
        renumber_pdbfile_to_aho(pdbfile,TAG,tmp_ahofile,True)
        inpose = rosetta.Pose()
        rosetta.core.import_pose.pose_from_pdb( inpose , tmp_ahofile )
        pdbseq = inpose.sequence()
        pos_cdr1s = inpose.pdb_info().pdb2pose(chainid,24)
        pos_cdr1e = inpose.pdb_info().pdb2pose(chainid,41)
        pos_cdr2s = inpose.pdb_info().pdb2pose(chainid,55)
        pos_cdr2e = inpose.pdb_info().pdb2pose(chainid,91)
        pos_cdr3s = inpose.pdb_info().pdb2pose(chainid,108)
        pos_cdr3e = inpose.pdb_info().pdb2pose(chainid,138)
        CDR1_match = False            
        if not ( (int(pos_cdr1s) == 0) or (int(pos_cdr1e) == 0) ):
            pdb_cdr1 = pdbseq[pos_cdr1s:pos_cdr1e-1]
            if (str(seq_cdr1) == str(pdb_cdr1)):
                CDR1_match=True
                with open(TAG+'_TCR_CDR1_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
                    fr.write(">"+record.id+"_"+TAG+"\n")
                    fr.write(pdb_cdr1+"\n")
        CDR2_match = False   
        if not ( (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) ):
            pdb_cdr2hv4 = pdbseq[pos_cdr2s:pos_cdr2e-1]
            if (str(seq_cdr2hv4) == str(pdb_cdr2hv4)):
                CDR2_match=True
                with open(TAG+'_TCR_CDR2HV4_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
                    fr.write(">"+record.id+"_"+TAG+"\n")
                    fr.write(pdb_cdr2hv4+"\n")
        CDR3_match = False
        if not ( (int(pos_cdr3s) == 0) or (int(pos_cdr3e) == 0) ):
            pdb_cdr3 = pdbseq[pos_cdr3s:pos_cdr3e-1]
            if (str(seq_cdr3) == str(pdb_cdr3)):
                CDR3_match=True
                with open(TAG+'_TCR_CDR3_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
                    fr.write(">"+record.id+"_"+TAG+"\n")
                    fr.write(pdb_cdr3+"\n")
        FW_match = False            
        if not ( (int(pos_cdr1s) == 0) or (int(pos_cdr1e) == 0) or (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) or (int(pos_cdr3s) == 0) or (int(pos_cdr3e) == 0) ):
            pdb_fw = pdbseq[pos_cdr1s-20:pos_cdr1s]+pdbseq[pos_cdr1e-1:pos_cdr2s]+pdbseq[pos_cdr2e-1:pos_cdr3s]+pdbseq[pos_cdr3e-1:pos_cdr3e+9]
            if (str(seq_fw) == str(pdb_fw)):
                FW_match=True
                with open(TAG+'_TCR_FW.fasta', 'a') as fr:
                    fr.write(">"+record.id+"_"+TAG+"\n")
                    fr.write(pdb_fw+"\n")
                    
        #create aho for pdb structure with missing residues
        if FW_match is True:
            pdbsubpose = graft.return_region( inpose, int(pos_cdr1s-19), int(pos_cdr3e+9) )
            pdbsubpose.dump_pdb(cap_pdbfile)
            graftpdbfname = change_chain_id(cap_pdbfile,"A",chainid)
            renumber_pdbfile_to_aho(graftpdbfname,TAG,cap_ahofile,True)
        else:
            if ( (CDR1_match is True) or (CDR2_match is True) or (CDR3_match is True) ):
                renumber_pdbfile_to_aho(pdbfile,TAG,nocap_ahofile,False)
                
