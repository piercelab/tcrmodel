#!/usr/bin/env python

import os
from glob import glob
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select, CaPPBuilder, PPBuilder
from Bio.PDB.StructureBuilder import StructureBuilder 
from Bio import SeqIO
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()

def seq_from_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id == inp_id:
            return rec.seq

def fasta_record_from_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id == inp_id:
            return rec

def fasta_record_from_trunc_id(inp_id, fasta_file):
    for rec in SeqIO.parse(fasta_file, "fasta"):
        if rec.id[:8] == inp_id[:8]:
            return rec

def get_unique_by_seq_and_resolution(curr_file):
    fasta_sequences = SeqIO.parse(open(curr_file),'fasta')
    ori_list = []
    for val in fasta_sequences:
        col = [val.id,val.seq,val.id[9:],val]
        ori_list.append(col)
    tmp_ori_list = ori_list
    unq_ori_list = []
    for item1 in ori_list:
        flag=True
        for unqitem in unq_ori_list:
            if ( str(item1[1]) == str(unqitem[1]) ):
                flag=False
        if flag:
            newitem = item1
            for item2 in tmp_ori_list:
                if ( str(newitem[1]) == str(item2[1]) ):
                    if (newitem[2] > item2[2]):
                        newitem = item2
            unq_ori_list.append(newitem)
    return unq_ori_list


def check_seq_match(inpseq,fasta_file):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for val in fasta_sequences:
        if str(inpseq.strip()) == str((val.seq).strip()):
            return True


def make_complex_from_chains(pdbchains,outfile=None):
    cwd = os.getcwd()
    pdb_dirpath = os.path.join(cwd,"../updated_templates/tcr/pdb")
    if not os.path.exists(pdb_dirpath):
        os.makedirs(pdb_dirpath)
    if outfile is None:
        #outfile = apdb[:4]+"_"+achainid+bchainid+"_aho.pdb"
        outfn = pdbchains[0][:4]+"_aho.pdb"
        outfile = os.path.join(pdb_dirpath,outfn)

    if len(pdbchains) < 2:
        #sb = parser.structure_builder.init_structure('pdb') 
        for pdbchain in pdbchains:
            chainid = pdbchain[5:6]
            pdbchainfile = pdbchain+"_aho.pdb"
            structure = parser.get_structure("chainid", pdbchainfile) 
        io.set_structure(structure)
        io.save(outfile)
    else:
        apdb = pdbchains[0]+"_aho.pdb"
        bpdb = pdbchains[1]+"_aho.pdb"
        achainid = apdb[5:6]
        astructure = parser.get_structure("Achain", apdb)
        bchainid = bpdb[5:6]
        bstructure = parser.get_structure("Bchain", bpdb)
        complexstruct = astructure
        bchain = bstructure[0][bchainid]
        newbchainid = bchainid
        if (achainid == bchainid):
            if achainid.isupper(): newbchainid = bchainid.lower()
            if achainid.islower(): newbchainid = bchainid.upper()
            bchain.id = newbchainid
        complexstruct[0].add(bchain)
        io.set_structure(complexstruct)
        io.save(outfile)
    return outfile

