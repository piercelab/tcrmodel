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


import sys
sys.path.insert(0,'/www/tcrmodel/update/tcr_template_structures/scripts')

import pyrosetta
pyrosetta.init(extra_options = "-mute basic -mute core -mute protocols -ignore_zero_occupancy false")
import pyrosetta.rosetta.protocols.grafting as graft
from pyrosetta import Pose
from pyrosetta import pose_from_file


from TCR_functions import *
script, fasta_file, TAG = argv

if TAG == 'A':
    aho_cdr1_s = 24
    aho_cdr1_e = 42
    aho_cdr2_s = 56
    aho_cdr2_e = 71
    aho_cdr2hv4_s = 56
    aho_cdr2hv4_e = 90
    aho_cdr3_s = 107
    aho_cdr3_e = 138
    aho_hv4_s = 81
    aho_hv4_e = 90

if TAG == 'B':
    aho_cdr1_s = 24
    aho_cdr1_e = 42
    aho_cdr2_s = 56
    aho_cdr2_e = 70
    aho_cdr2hv4_s = 56
    aho_cdr2hv4_e = 90
    aho_cdr3_s = 107
    aho_cdr3_e = 138
    aho_hv4_s = 81
    aho_hv4_e = 90

cap_begin = 20
cap_end = 10

def extract_region(ahopose, pdb_start_pos, pdb_end_pos, chainid, outfile):
    pose_start_pos = ahopose.pdb_info().pdb2pose(chainid, pdb_start_pos)
    pose_end_pos = ahopose.pdb_info().pdb2pose(chainid, pdb_end_pos)
    #print pdb_start_pos, pdb_end_pos, pose_start_pos, pose_end_pos                                                                     
    if (pose_start_pos == 0) or (pose_end_pos == 0):
        return "ERROR"
    else:
        cdrpose = Pose()
        cdrpose = graft.return_region( ahopose, pose_start_pos, pose_end_pos)
        cdrpose.dump_pdb(outfile)
        return cdrpose

cdrfilename = TAG+"_CDRDATA.txt"
if os.path.exists(cdrfilename):
    os.remove(cdrfilename)
cdrfile = open(cdrfilename, "a")
if TAG == "A": domtag = "Alpha"
if TAG == "B": domtag = "Beta"
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    pdbid = record.id[:4]
    chainid = record.id[5:6]
    resolution = extract_resolution_from_pdb(record.id)
    outtag = pdbid+"_"+chainid+"_"+TAG
    pdbfile = outtag+".pdb"
    cap_pdbfile = outtag+"_cap.pdb"
    cap_ahofile = outtag+"_aho.pdb"
    nocap_ahofile = outtag+"_aho.pdb"
    outheader = pdbid+"_"+chainid+"_"+TAG+"_"+str(resolution)
    print outtag, outheader
    #print record.seq
    #seq based
    fseq = record.seq
    if check_tcr_alpha_beta_domain(record.seq) is True:
        anarci_fseq = get_tcr_domain_seq(record.seq,TAG)
    else:
        print pdbid, chainid, "WARNING: tcr_alpha_beta_domain not found with ANARCI"
        
    cdrfromseq = get_cdr_from_seq_by_aho_num(record.seq,TAG)
    seq_cdr1_aho = str(cdrfromseq[0])
    seq_cdr2_aho = str(cdrfromseq[1])  
    seq_cdr2hv4_aho = str(cdrfromseq[2])
    seq_hv4_aho = str(cdrfromseq[3])  
    seq_cdr3_aho = str(cdrfromseq[4])  
    
    regexres = assign_CDRs_using_REGEX(anarci_fseq, TAG)
    if regexres:
        seq_vdomain = regexres.groups()[0]
        seq_gm = regexres.groups()[1] + regexres.groups()[2] + regexres.groups()[3] + regexres.groups()[4] + regexres.groups()[5] + regexres.groups()[7]
        seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
        seq_cdr1 = regexres.groups()[2]
        seq_cdr2hv4 = regexres.groups()[4]
        seq_cdr3 = regexres.groups()[6]   
        #print "cap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
        with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_vdomain+"\n")
    else:
        regexres = assign_CDRs_using_REGEX(fseq, TAG)
        if regexres:
            seq_vdomain = regexres.groups()[0]
            seq_gm = regexres.groups()[1] + regexres.groups()[2] + regexres.groups()[3] + regexres.groups()[4] + regexres.groups()[5] + regexres.groups()[7]
            seq_fw = regexres.groups()[1] + regexres.groups()[3] + regexres.groups()[5] + regexres.groups()[7]
            seq_cdr1 = regexres.groups()[2]
            seq_cdr2hv4 = regexres.groups()[4]
            seq_cdr3 = regexres.groups()[6]   
            #print "cap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
            with open(TAG+'_TCR_VDOMAIN.fasta', 'a') as fr:
                fr.write(">"+outheader+"\n")
                fr.write(seq_vdomain+"\n")
            
        else:
            regexres = assign_CDRs_using_REGEX(fseq, TAG, True)
            if regexres:
                seq_vdomain = "None"
                seq_fw = "None"
                seq_gm = "None"
                seq_cdr1 = regexres.groups()[2]
                seq_cdr2hv4 = regexres.groups()[4]
                seq_cdr3 = regexres.groups()[6]
            #print "nocap", seq_vdomain, seq_fw, seq_cdr1, seq_cdr2hv4, seq_cdr3
            else:
                print outtag, "FAIL: RegEx failed to identify TCR domain"
                continue

    pdbgzfile = os.path.join("/www/tcrmodel/update/tcr_template_structures/rcsb_update/structures/"+pdbid+".pdb.gz")
    get_pdb_from_gzpdbfile(record.id[:4],record.id[5:6],pdbfile,pdbgzfile)
    strseq = pdb_to_fasta(pdbfile)
    pdbregex = re.search(seq_vdomain,strseq)
    if pdbregex:
        capbegin = pdbregex.span()[0]
        capend = pdbregex.span()[1]
        #print capbegin, capend
        cap_pdb(pdbfile,capbegin,capend,cap_pdbfile)
        renumber_pdbfile_to_aho(cap_pdbfile,TAG,cap_ahofile,True)
        with open(TAG+'_TCR_VDOMAIN_cap.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_vdomain+"\n")
        with open(TAG+'_TCR_GM.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_gm+"\n")
        with open(TAG+'_TCR_FW.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_fw+"\n")
        with open(TAG+'_TCR_CDR1_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr1+"\n")
        with open(TAG+'_TCR_CDR2_'+str(len(seq_cdr2_aho))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr2_aho+"\n")
        with open(TAG+'_TCR_CDR2HV4_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr2hv4+"\n")
        with open(TAG+'_TCR_HV4_'+str(len(seq_hv4_aho))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_hv4_aho+"\n")
        with open(TAG+'_TCR_CDR3_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
            fr.write(">"+outheader+"\n")
            fr.write(seq_cdr3+"\n")

        '''
        capahopose = Pose()
        capahopose = pose_from_file( cap_ahofile )

        outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr1"+"_"+str(len(seq_cdr1))+".pdb"
        extract_region(capahopose, aho_cdr1_s, aho_cdr1_e, chainid, outfile)
        cdrfile.write("cdr1"+" "+seq_cdr1+" "+pdbid+" "+chainid+" "+str(len(seq_cdr1))+" "+domtag+" "+"OK"+"\n")

        outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr2"+"_"+str(len(seq_cdr2_aho))+".pdb"
        extract_region(capahopose, aho_cdr2_s, aho_cdr2_e, chainid, outfile)
        cdrfile.write("cdr2"+" "+seq_cdr2_aho+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2_aho))+" "+domtag+" "+"OK"+"\n")

        outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr2hv4"+"_"+str(len(seq_cdr2hv4))+".pdb"
        extract_region(capahopose, aho_cdr2hv4_s, aho_cdr2hv4_e, chainid, outfile)
        cdrfile.write("cdr2hv4"+" "+seq_cdr2hv4+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2hv4))+" "+domtag+" "+"OK"+"\n")

        outfile = pdbid+"_"+chainid+"_"+TAG+"_hv4"+"_"+str(len(seq_hv4_aho))+".pdb"
        extract_region(capahopose, aho_hv4_s, aho_hv4_e, chainid, outfile)
        cdrfile.write("hv4"+" "+seq_hv4_aho+" "+pdbid+" "+chainid+" "+str(len(seq_hv4_aho))+" "+domtag+" "+"OK"+"\n")

        outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr3"+"_"+str(len(seq_cdr3))+".pdb"
        extract_region(capahopose, aho_cdr3_s, aho_cdr3_e, chainid, outfile)
        cdrfile.write("cdr3"+" "+seq_cdr3+" "+pdbid+" "+chainid+" "+str(len(seq_cdr3))+" "+domtag+" "+"OK"+"\n")
        '''
    else:
        tmp_ahofile = "temp.aho.pdb"
        renumber_pdbfile_to_aho(pdbfile,TAG,tmp_ahofile,True)
        inpose = Pose()
        inpose = pose_from_file( tmp_ahofile )
        pdbseq = inpose.sequence()

        pdb_cdr1=""
        pdb_cdr2=""
        pdb_cdr2hv4=""
        pdb_hv4=""
        pdb_cdr3=""
        pos_cdr1s = inpose.pdb_info().pdb2pose(chainid,aho_cdr1_s)
        pos_cdr1e = inpose.pdb_info().pdb2pose(chainid,aho_cdr1_e)
        pos_cdr2s = inpose.pdb_info().pdb2pose(chainid,aho_cdr2_s)
        pos_cdr2e = inpose.pdb_info().pdb2pose(chainid,aho_cdr2_e)
        pos_cdr2hv4s = inpose.pdb_info().pdb2pose(chainid,aho_cdr2hv4_s)
        pos_cdr2hv4e = inpose.pdb_info().pdb2pose(chainid,aho_cdr2hv4_e)
        pos_cdr3s = inpose.pdb_info().pdb2pose(chainid,aho_cdr3_s)
        pos_cdr3e = inpose.pdb_info().pdb2pose(chainid,aho_cdr3_e)
        pos_hv4s = inpose.pdb_info().pdb2pose(chainid,aho_hv4_s)
        pos_hv4e = inpose.pdb_info().pdb2pose(chainid,aho_hv4_e)

        CDR1_match = False
        if not ( (int(pos_cdr1s) == 0) or (int(pos_cdr1e) == 0) ):
            pdb_cdr1 = pdbseq[pos_cdr1s-1:pos_cdr1e]
            if (str(seq_cdr1) == str(pdb_cdr1)):
                CDR1_match=True
                with open(TAG+'_TCR_CDR1_noclean_'+str(len(seq_cdr1))+'.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_cdr1+"\n")
                outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr1"+"_"+str(len(seq_cdr1))+".pdb"
                extract_region(inpose, aho_cdr1_s, aho_cdr1_e, chainid, outfile)
                cdrfile.write("cdr1"+" "+seq_cdr1+" "+pdbid+" "+chainid+" "+str(len(seq_cdr1))+" "+domtag+" "+"OK"+"\n")
            else:
                cdrfile.write("cdr1"+" "+seq_cdr1+" "+pdbid+" "+chainid+" "+str(len(seq_cdr1))+" "+domtag+" "+"Mismatch"+"\n")
        else:
            cdrfile.write("cdr1"+" "+seq_cdr1+" "+pdbid+" "+chainid+" "+str(len(seq_cdr1))+" "+domtag+" "+"Mismatch"+"\n")


        CDR2_match = False
        if not ( (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) ):
            pdb_cdr2 = pdbseq[pos_cdr2s-1:pos_cdr2e]
            if (str(seq_cdr2_aho) == str(pdb_cdr2)):
                CDR1_match=True
                with open(TAG+'_TCR_CDR2_noclean_'+str(len(seq_cdr2_aho))+'.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_cdr1+"\n")
                outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr2"+"_"+str(len(seq_cdr2_aho))+".pdb"
                extract_region(inpose, aho_cdr2_s, aho_cdr2_e, chainid, outfile)
                cdrfile.write("cdr2"+" "+seq_cdr2_aho+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2_aho))+" "+domtag+" "+"OK"+"\n")
            else:
                cdrfile.write("cdr2"+" "+seq_cdr2_aho+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2_aho))+" "+domtag+" "+"Mismatch"+"\n")
        else:
            cdrfile.write("cdr2"+" "+seq_cdr2_aho+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2_aho))+" "+domtag+" "+"Mismatch"+"\n")

        CDR2HV4_match = False   
        if not ( (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) ):
            pdb_cdr2hv4 = pdbseq[pos_cdr2hv4s-1:pos_cdr2hv4e]
            if (str(seq_cdr2hv4) == str(pdb_cdr2hv4)):
                CDR2HV4_match=True
                with open(TAG+'_TCR_CDR2HV4_noclean_'+str(len(seq_cdr2hv4))+'.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_cdr2hv4+"\n")
                outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr2hv4"+"_"+str(len(seq_cdr2hv4))+".pdb"
                extract_region(inpose, aho_cdr2hv4_s, aho_cdr2hv4_e, chainid, outfile)
                cdrfile.write("cdr2hv4"+" "+seq_cdr2hv4+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2hv4))+" "+domtag+" "+"OK"+"\n")
            else:
                cdrfile.write("cdr2hv4"+" "+seq_cdr2hv4+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2hv4))+" "+domtag+" "+"Mismatch"+"\n")
        else:
            cdrfile.write("cdr2hv4"+" "+seq_cdr2hv4+" "+pdbid+" "+chainid+" "+str(len(seq_cdr2hv4))+" "+domtag+" "+"Mismatch"+"\n")

        HV4_match = False
        if not ( (int(pos_hv4s) == 0) or (int(pos_hv4e) == 0) ):
            pdb_hv4 = pdbseq[pos_hv4s-1:pos_hv4e]
            if (str(seq_hv4_aho) == str(pdb_hv4)):
                HV4_match=True
                with open(TAG+'_TCR_HV4_noclean_'+str(len(seq_hv4_aho))+'.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_hv4+"\n")
                outfile = pdbid+"_"+chainid+"_"+TAG+"_hv4"+"_"+str(len(seq_hv4_aho))+".pdb"
                extract_region(inpose, aho_hv4_s, aho_hv4_e, chainid, outfile)
                cdrfile.write("hv4"+" "+seq_hv4_aho+" "+pdbid+" "+chainid+" "+str(len(seq_hv4_aho))+" "+domtag+" "+"OK"+"\n")
            else:
                cdrfile.write("hv4"+" "+seq_hv4_aho+" "+pdbid+" "+chainid+" "+str(len(seq_hv4_aho))+" "+domtag+" "+"Mismatch"+"\n")
        else:
            cdrfile.write("hv4"+" "+seq_hv4_aho+" "+pdbid+" "+chainid+" "+str(len(seq_hv4_aho))+" "+domtag+" "+"Mismatch"+"\n")

        CDR3_match = False
        if not ( (int(pos_cdr3s) == 0) or (int(pos_cdr3e) == 0) ):
            pdb_cdr3 = pdbseq[pos_cdr3s-1:pos_cdr3e]
            if (str(seq_cdr3) == str(pdb_cdr3)):
                CDR3_match=True
                with open(TAG+'_TCR_CDR3_noclean_'+str(len(seq_cdr3))+'.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_cdr3+"\n")
                outfile = pdbid+"_"+chainid+"_"+TAG+"_cdr3"+"_"+str(len(seq_cdr3))+".pdb"
                extract_region(inpose, aho_cdr3_s, aho_cdr3_e, chainid, outfile)
                cdrfile.write("cdr3"+" "+seq_cdr3+" "+pdbid+" "+chainid+" "+str(len(seq_cdr3))+" "+domtag+" "+"OK"+"\n")
            else:
                cdrfile.write("cdr3"+" "+seq_cdr3+" "+pdbid+" "+chainid+" "+str(len(seq_cdr3))+" "+domtag+" "+"Mismatch"+"\n")
                print "cdr3", "Mismatch", pdb_cdr3, seq_cdr3
        else:
            cdrfile.write("cdr3"+" "+seq_cdr3+" "+pdbid+" "+chainid+" "+str(len(seq_cdr3))+" "+domtag+" "+"Mismatch"+"\n")
            print "cdr3", "Mismatch", pdb_cdr3, seq_cdr3            

        FW_match = False            
        if not ( (int(pos_cdr1s) == 0) or (int(pos_cdr1e) == 0) or (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) or (int(pos_cdr3s) == 0) or (int(pos_cdr3e) == 0) ):
            pdb_fw = pdbseq[pos_cdr1s-1-cap_begin:pos_cdr1s-1]+pdbseq[pos_cdr1e:pos_cdr2hv4s-1]+pdbseq[pos_cdr2hv4e:pos_cdr3s-1]+pdbseq[pos_cdr3e:pos_cdr3e+cap_end]
            if (str(seq_fw) == str(pdb_fw)):
                FW_match=True
                with open(TAG+'_TCR_FW_noclean.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_fw+"\n")
            #print "TEST fw : ", seq_fw, pdb_fw

        GM_match = False            
        if not ( (int(pos_cdr1s) == 0) or (int(pos_cdr1e) == 0) or (int(pos_cdr2s) == 0) or (int(pos_cdr2e) == 0) or (int(pos_cdr3s) == 0) or (int(pos_cdr3e) == 0) ):
            pdb_gm = pdbseq[pos_cdr1s-1-cap_begin:pos_cdr3s-1]+pdbseq[pos_cdr3e:pos_cdr3e+cap_end]
            if (str(seq_gm) == str(pdb_gm)):
                GM_match=True
                with open(TAG+'_TCR_GM_noclean.fasta', 'a') as fr:
                    fr.write(">"+outheader+"\n")
                    fr.write(pdb_gm+"\n")

        #create aho for pdb structure with missing residues
        if FW_match is True:
            pdbsubpose = graft.return_region( inpose, int(pos_cdr1s-cap_begin), int(pos_cdr3e+cap_end) )
            pdbsubpose.dump_pdb(cap_pdbfile)
            graftpdbfname = change_chain_id(cap_pdbfile,"A",chainid)
            renumber_pdbfile_to_aho(graftpdbfname,TAG,cap_ahofile,True)
        else:
            if ( (CDR1_match is True) or (CDR2HV4_match is True) or (CDR3_match is True) ):
                renumber_pdbfile_to_aho(pdbfile,TAG,nocap_ahofile,False)
                
cdrfile.close()
