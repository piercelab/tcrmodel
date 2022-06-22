#!/usr/bin/python 

import os, sys
import glob
import subprocess
from Bio import SeqIO
from subprocess import Popen, PIPE

sys.path.append('/home/rg/scripts')
from TCRmodeller_functions import *


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
            rmsd = line.split(' ')[-1].rstrip()
            return rmsd


def find_best_template(inpcdr, all_cdr_fasta_files):
    best_score = -99999
    best_seq = ""
    best_seqid = ""
    for cdrfile in all_cdr_fasta_files:
        for record in SeqIO.parse(cdrfile, "fasta"):
            if not len(inpcdr.seq) == len(record.seq):
                continue
            if inpcdr.seq == record.seq:
                continue
            print calculate_identity_score(inpcdr.seq, record.seq)
            if calculate_identity_score(inpcdr.seq, record.seq) > 90:
                continue
            score = score_alignment(inpcdr.seq, record.seq)
            if (score > best_score):
                best_score = score
                best_seq = record.seq
                best_seqid = record.id
                
    cdr_pdb_alpha_dir = "/home/rg/cdr_alpha"
    cdr_pdb_beta_dir = "/home/rg/cdr_beta"
    
    reference = os.path.join(cdr_pdb_alpha_dir, inpcdr.id + ".pdb")
    if not os.path.isfile(reference):
        reference = os.path.join(cdr_pdb_beta_dir, inpcdr.id + ".pdb")
    
    mobile = os.path.join(cdr_pdb_alpha_dir, best_seqid + ".pdb")
    if not os.path.isfile(mobile):
        mobile = os.path.join(cdr_pdb_beta_dir, best_seqid + ".pdb")

    rms = get_rms("profit.in", reference, mobile)
    #print inpcdr.seq, best_seq, inpcdr.id, best_seqid, best_score, rms
    if not best_score == -99999:
        print str(best_score)+"\t"+str(rms)



def find_best_template_by_percent_identity(inpcdr, all_cdr_fasta_files):
    best_score = -99999
    best_seq = ""
    best_seqid = ""
    for cdrfile in all_cdr_fasta_files:
        for record in SeqIO.parse(cdrfile, "fasta"):
            if inpcdr.seq == record.seq:
                continue
            score = calculate_identity_score(inpcdr.seq, record.seq)
            if (score > best_score):
                best_score = score
                best_seq = record.seq
                best_seqid = record.id
                
    cdr_pdb_alpha_dir = "/home/rg/cdr_alpha"
    cdr_pdb_beta_dir = "/home/rg/cdr_beta"
    
    reference = os.path.join(cdr_pdb_alpha_dir, inpcdr.id + ".pdb")
    if not os.path.isfile(reference):
        reference = os.path.join(cdr_pdb_beta_dir, inpcdr.id + ".pdb")
    
    mobile = os.path.join(cdr_pdb_alpha_dir, best_seqid + ".pdb")
    if not os.path.isfile(mobile):
        mobile = os.path.join(cdr_pdb_beta_dir, best_seqid + ".pdb")

    rms = get_rms("profit.in", reference, mobile)
    print inpcdr.seq, best_seq, inpcdr.id, best_seqid, best_score, rms
    if not best_score == -99999:
        print str(best_score)+"\t"+str(rms)





# main

from sys import argv
script, cdr_fasta_file = argv

myPath= os.getcwd()

profit_infile = 'profit.in'
create_profit_infile(profit_infile)
all_cdr_fasta_files = glob.glob("/home/rg/cdr_template_coverage/clear_*.fa")

fasta_sequences = SeqIO.parse(open(cdr_fasta_file),'fasta')
for inpcdr in fasta_sequences:
    find_best_template(inpcdr, all_cdr_fasta_files)
    #find_best_template_by_percent_identity(inpcdr, all_cdr_fasta_files)

print "done"

#print score_alignment("AFMDSNYQLI", "ASMDSNYQLI")
#print score_alignment("AFMDSNYQLII","ASMDSNYQLII")
