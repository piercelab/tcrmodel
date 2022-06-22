#!/usr/bin/env python

import sys, os, re
import gzip
from sys import argv
from Bio.PDB import PDBParser, Selection, PDBIO, CaPPBuilder, Select
from Bio import SeqIO

parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()

import rosetta
rosetta.init(extra_options = "-database /TCRmodeller/programs/PyRosetta/rosetta_database -mute basic -mute core -mute protocols -renumber_pdb -ignore_unrecognized_res -ignore_zero_occupancy false")
import rosetta.protocols.grafting as graft

script, pdbcode, chainida, chainidb = argv

def find_tcr_vdomain_seq(seq, tag):
    if tag == "alpha":
        regex = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[YF][A-Z]{13}([A-Z]{6,11})[A-Z]{15,30}[DL][A-Z]{2,3}Y[A-Z][CW][A-Z]([A-Z]{7,16}[FW])G[A-Z]G[A-Z]{0,7}[PA]*"
    if tag == "beta":
        regex = "[A-Z]{0,23}C[A-Z]([A-Z]{8,12}W)[Y][A-Z]{13}([A-Z]{6,11})[A-Z]{15,40}[YLF][A-Z][CW][A-Z]([A-Z]{7,17}[F])G[A-Z]G[A-Z]{0,7}[E]*"
    res = re.search(regex, str(seq))
    if res:
        return res
        #print res.group(), res.start(), res.end()
    else:
        print None


class nonHetSelect(Select):
    def accept_residue(self,residue):
        if residue.id[0] == ' ':
            return 1
        else:
            return 0

from TCR_functions import *
gzpdbfile_path =  "/TCRmodeller/PDB_RELEASE/pdb_structures" +  '/%s/pdb%s.ent.gz' %(pdbcode[1:3], pdbcode) 
#gzpdbfile_path =  pdbcode+".pdb.gz" 
#import wget
#url = "http://www.rcsb.org/pdb/files/"+gzpdbfile_path
#filename = wget.download(url)
gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
pdbfile = parser.get_structure("PDB", gzpdbfile)

mychaina = pdbfile[0][chainida]
io.set_structure(mychaina)
io.save('tmpa.pdb', nonHetSelect())
achainseq = ""
for ppe in ppb.build_peptides(mychaina):
    achainseq += str(ppe.get_sequence())

mychainb = pdbfile[0][chainidb]
io.set_structure(mychainb)
io.save('tmpb.pdb', nonHetSelect())
bchainseq = ""
for ppe in ppb.build_peptides(mychainb):
    bchainseq += str(ppe.get_sequence())

fasta_file =  "/TCRmodeller/PDB_RELEASE/pdb_sequences/pdb_seqres.txt"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:
    if record.id[:4] == pdbcode and record.id[5:6] == chainida:
        afasta = record.seq 
    if record.id[:4] == pdbcode and record.id[5:6] == chainidb:
        bfasta = record.seq

#print achainseq
#print afasta

#print bchainseq
#print bfasta


print "\ntcrdomain"

print achainseq
print afasta
#tmp1a = find_variable_domain_using_regex(achainseq, 'a')
tmp1a = assign_CDRs_using_REGEX(achainseq, 'A')
tmp1agr = tmp1a.groups()
#tmp2a = find_variable_domain_using_regex(afasta, 'a')
tmp2a = assign_CDRs_using_REGEX(afasta, 'A')
tmp2agr = tmp2a.groups()
 
print tmp1agr
print tmp2agr   

#tmp1b = find_variable_domain_using_regex(bchainseq, 'b')
tmp1b = assign_CDRs_using_REGEX(bchainseq, 'B')
tmp1bgr = tmp1b.groups()
tmp2b = assign_CDRs_using_REGEX(bchainseq, 'B')
#tmp2b = find_variable_domain_using_regex(bfasta, 'b')
tmp2bgr = tmp2b.groups()

if re.search(tmp1agr[0], tmp1agr[0]) and re.search(tmp1bgr[0], tmp2bgr[0]):
    print "RES MATCH", pdbcode
else:
    print "RES NOMATCH", pdbcode


'''
astrres = assign_CDRs_using_REGEX(achainseq, 'A')
atcr = astrres.groups()
aseqres = assign_CDRs_using_REGEX(afasta, 'A')
aseqtcr = aseqres.groups()
if atcr[0] == aseqtcr[0]:
    print "A Match", atcr[0] 
    print "A Match", aseqtcr[0]
else:
    print "A NoMatch", atcr[0]
    print "A NoMatch", aseqtcr[0]

bstrres = assign_CDRs_using_REGEX(bchainseq, 'B')
btcr = bstrres.groups()
bseqres = assign_CDRs_using_REGEX(bfasta, 'B')
bseqtcr = bseqres.groups()
if btcr[0] == bseqtcr[0]:
    print "B Match", btcr[0] 
    print "B Match", bseqtcr[0]
else:
    print "B NoMatch", btcr[0]
    print "B NoMatch", bseqtcr[0]


atruncdomain = atcr[0]
fra = atcr[1]+atcr[3]+atcr[5]+atcr[7]
cdr1a = atcr[2]
cdr2a = atcr[4]
cdr3a = atcr[6]
atcrs = [cdr1a,cdr2a,cdr3a,atruncdomain,fra]
cdr1a_span = ares.span(3)
cdr2a_span = ares.span(5)
cdr3a_span = ares.span(7)

bres = assign_CDRs_using_REGEX(bchainseq, 'B')
btcr = bres.groups()
btruncdomain = btcr[0]
frb = btcr[1]+btcr[3]+btcr[5]+btcr[7]
cdr1b = btcr[2]
cdr2b = btcr[4]
cdr3b = btcr[6]
btcrs = [cdr1b,cdr2b,cdr3b,btruncdomain,frb]
cdr1b_span = bres.span(3)
cdr2b_span = bres.span(5)
cdr3b_span = bres.span(7)


        
protposea = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( protposea , 'tmpa.pdb' )
apose = rosetta.Pose()
apose = graft.return_region( protposea, resa.start()+1, resa.end())
apose.dump_pdb("apose.pdb")    

protposeb = rosetta.Pose()
rosetta.core.import_pose.pose_from_pdb( protposeb , 'tmpb.pdb' )
apose = rosetta.Pose()
apose = graft.return_region( protposeb, resb.start()+1, resb.end())
apose.dump_pdb("bpose.pdb")    

outfilename = pdbcode + '.trunc.pdb'
OUT = open(outfilename, 'w+')

INA = open("apose.pdb")
for line in INA.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
            newline  = line[:21] + "A" + line[22:]
            OUT.write(newline)
OUT.write('TER\n')
INA.close()

INB = open("bpose.pdb")
for line in INB.readlines():
    if line[0:4] == 'ATOM':
        if line[21:22] == 'A':
            newline  = line[:21] + "B" + line[22:]
            OUT.write(newline)
OUT.write('TER\n')
INB.close()
OUT.close()



seqresa = find_tcr_vdomain_seq(afasta, "alpha")
print seqresa.group(), seqresa.start(), seqresa.end()
atcr = seqresa.groups()
atruncdomain = atcr[0]
cdr1a = atcr[1]
cdr2a = atcr[2]
cdr3a = atcr[3]
atcrs = [cdr1a,cdr2a,cdr3a,atruncdomain,fra]
cdr1a_span = ares.span(3)
cdr2a_span = ares.span(5)
cdr3a_span = ares.span(7)

seqresb = find_tcr_vdomain_seq(bfasta, "beta")
print seqresb.group(), seqresb.start(), seqresb.end()

'''
