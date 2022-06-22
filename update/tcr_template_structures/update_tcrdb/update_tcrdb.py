#!/usr/bin/python

import os
import sys
import subprocess 
import sqlite3
import json
import urllib2

from Bio import SeqIO
sys.path.insert(0,'/home/gowthamanr/tcr_update/template_update_scripts')
from TCR_functions import *


#global
trav_regex = [
    "^[A-Z]*([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y][A-Z])[C][A-Z]*",
    "^[A-Z]*([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y][A-Z])[W][A-Z]*"
    #"^[A-Z]*([C][A-Z]{1,19}[W][A-Z]{1,70}[Y][A-Z])[C][A-Z]*", for shorter seqs
]

trbv_regex = [
    "([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y|F|L][A-Z])[C][A-Z]*",
    "^[A-Z]*([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y|F][A-Z])CTC[A-Z]*"
]

cwd = os.getcwd()
os.chdir("../rcsb_update") 
#subprocess.call(["./update.sh"])
tcr_alpha_chain = "../rcsb_update/tcr_alpha.fasta"
tcr_beta_chain = "../rcsb_update/tcr_beta.fasta"
#tcr_beta_chain = "tmp.fasta"
all_pdb_seq = "../rcsb_update/pdb_seqres.txt"
os.chdir(cwd) 

def update_complexes():
    ab_chain_list = []
    for record in SeqIO.parse(tcr_alpha_chain, "fasta"):
        pdbid = record.id[:4]
        if pdbid not in ab_chain_list:
            ab_chain_list.append(pdbid)
    for record in SeqIO.parse(tcr_beta_chain, "fasta"):
        pdbid = record.id[:4]
        if pdbid not in ab_chain_list:
            ab_chain_list.append(pdbid)

    all_ab_seq = "fasta_files/all_ab_seq.fasta"
    with open(all_ab_seq, "w") as output_handle:
        for record in SeqIO.parse(all_pdb_seq, "fasta"):
            pdbid = record.id[:4]
            if pdbid in ab_chain_list:
                SeqIO.write(record, output_handle, "fasta")


def calculate_identity_score( query, content):
    if (len(query) != len(content)): return -99999
    percent_identity = 0
    num_simi_res = 0
    for x in xrange(0, len(query)):
        if ( query[x] == content[x] ): num_simi_res += 1
    percent_identity = ( num_simi_res / float(len(query)) ) * 100;
    return percent_identity;

def get_truncated_trv_regex(inseq, regex_list):
    truncated_seq = None
    for tr_regex in regex_list:
        res = re.search(tr_regex, str(inseq))
        if res:
            truncated_seq = res.groups()[0]
            break
    return truncated_seq

def update_germline_gene(tcr_chain, TAG, regex_list, human_json, mouse_json):
    for record in SeqIO.parse(tcr_chain, "fasta"):
        pdbid_val=record.id[0:4]
        chain_val=record.id[5:6]
        print "adding ", pdbid_val, chain_val 
        cur.execute( """SELECT pdbid, chain FROM germline_gene WHERE pdbid=? AND chain=? AND tcr_type=?""",(pdbid_val, chain_val, TAG) )
        result = cur.fetchone()
        if result:
            # Record already exists
            continue
        else:
            print "adding ", pdbid_val, chain_val 
            curr_seq = get_truncated_trv(record.seq)
            print curr_seq
            curr_seq1 = get_truncated_trv_regex(record.seq, regex_list)
            print curr_seq1
            if curr_seq:
                my_subgroup_name = ''
                my_gene_name = ''
                my_allele_name = ''
                my_species = ''
                my_exact_match = ''
                match = False
                with open(human_json) as json_data:
                    hd = json.load(json_data)
                    for line in hd:
                        if ( str(line['seq']) == str(curr_seq) ):
                            match = True
                            my_subgroup_name = line['subgroup_name']
                            my_gene_name = line['gene_name']
                            my_allele_name = line['allele_name']
                            my_species = "Human"
                            my_exact_match = "MATCH"
                            break
                if not match:
                    with open(mouse_json) as json_data:
                        md = json.load(json_data)
                        for line in md:
                            if ( str(line['seq']) == str(curr_seq) ):
                                match = True
                                my_subgroup_name = line['subgroup_name']
                                my_gene_name = line['gene_name']
                                my_allele_name = line['allele_name']
                                my_species = "Mouse"
                                my_exact_match = "MATCH"
                                break
                if not match:                
                    best_score = -99999
                    best_subgroup_name = ''
                    best_gene_name = ''
                    best_allele_name = ''
                    species = ''
                    for line in hd:
                        score = calculate_identity_score(str(curr_seq), str(line['seq']))
                        if score > best_score:
                            best_score = score
                            best_subgroup_name = line['subgroup_name']
                            best_gene_name = line['gene_name']
                            best_allele_name = line['allele_name']
                            species = "Human"
                    for line in md:
                        score = calculate_identity_score(str(curr_seq), str(line['seq']))
                        if score > best_score:
                            best_score = score
                            best_subgroup_name = line['subgroup_name']
                            best_gene_name = line['gene_name']
                            best_allele_name = line['allele_name']
                            species = "Mouse"
                    if best_score > -99999:
                        my_subgroup_name = best_subgroup_name
                        my_gene_name = best_gene_name
                        my_allele_name = best_allele_name
                        my_species = species
                        my_exact_match = "NO_MATCH"
                    else:
                        print str(curr_seq)
                        print "germline_gene table : No germline match found , skipping pdb ", record.id
                        my_subgroup_name = "NA"
                        my_gene_name = "NA"
                        my_allele_name = "NA"
                        my_species = "NA"
                        my_exact_match = "NA"
            else:
                my_subgroup_name = "NA"
                my_gene_name = "NA"
                my_allele_name = "NA"
                my_species = "NA"
                my_exact_match = "NA"
                print "germline_gene table : failed to identify truncated domain, skipping pdb ", record.id
        cur.execute('''INSERT INTO germline_gene(pdbid,chain,tcr_type,subgroup_name,gene_name,allele_name,species,exact_match) VALUES (?, ?, ?, ?, ?, ?, ?, ?)''', (pdbid_val,chain_val,TAG,my_subgroup_name,my_gene_name,my_allele_name,my_species,my_exact_match) )
    return


def get_pdb_release_date(pdbid):
    xmlfile = urllib2.urlopen('https://www.rcsb.org/pdb/rest/customReport.xml?pdbids='+pdbid+'&customReportColumns=releaseDate')
    xmldata = xmlfile.read()
    xmlfile.close()
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xmldata)
    release_date = root[0][1].text
    return release_date

def update_tcr_chains(tcr_chain, TAG, longtag ):
    for record in SeqIO.parse(tcr_chain, "fasta"):
        pdbid_val=record.id[0:4]
        chain_val=record.id[5:6]
        cur.execute( """SELECT pdbid, chain FROM tcr_chains WHERE pdbid=? AND chain=?""",(pdbid_val, chain_val) )
        result = cur.fetchone()
        if result:
            # Record already exists
            continue
        print "adding ", pdbid_val, chain_val 
        desclist = []
        lenlist = []
        seq_val=str(record.seq)
        tcr_type_val=longtag
        description=record.description
        desclist = description.split(" ",3)
        desc_val=desclist[3]
        lenlist = str(desclist[2]).split(":")
        len_val = int(lenlist[1])
        release_date = get_pdb_release_date(pdbid_val)
        cur.execute( "INSERT OR IGNORE INTO tcr_chains(pdbid, chain, length, seq, tcr_type, desc, release_date) VALUES(?,?,?,?,?,?,?)",(pdbid_val, chain_val, len_val, seq_val, tcr_type_val, desc_val, release_date) )

def update_cdr_data(tcr_chain, TAG, longtag ):
    for record in SeqIO.parse(tcr_chain, "fasta"):
        pdbid_val=record.id[0:4]
        chain_val=record.id[5:6]
        cur.execute( """SELECT pdbid, chain FROM cdr_data WHERE pdbid=? AND chain=?""",(pdbid_val, chain_val) )
        result = cur.fetchone()
        if result:
            # Record already exists
            continue
        else:
            curr_seq = get_truncated_trv(record.seq)
            print "adding ", pdbid_val, chain_val 
            tcr_type = longtag
            remarks = ""
            cdr_from_seq_anarci = get_cdr_from_seq_by_aho_num_ext(record.seq,TAG)
            cdr_type = "cdr1"
            cdrseq = cdr_from_seq_anarci[0]
            cdr_len = len(str(cdrseq))
            cur.execute("INSERT INTO cdr_data (cdr_type, cdrseq, pdbid, chain, cdr_length, tcr_type, remarks) VALUES (?,?,?,?,?,?,?)", (cdr_type, cdrseq, pdbid_val, chain_val, cdr_len, tcr_type, remarks))
            
            cdr_type = "cdr2"
            cdrseq = cdr_from_seq_anarci[1]
            cdr_len = len(str(cdrseq))
            cur.execute("INSERT INTO cdr_data (cdr_type, cdrseq, pdbid, chain, cdr_length, tcr_type, remarks) VALUES (?,?,?,?,?,?,?)", (cdr_type, cdrseq, pdbid_val, chain_val, cdr_len, tcr_type, remarks))
            
            cdr_type = "cdr3_extnd"
            cdrseq = cdr_from_seq_anarci[4]
            cdr_len = len(str(cdrseq))
            cur.execute("INSERT INTO cdr_data (cdr_type, cdrseq, pdbid, chain, cdr_length, tcr_type, remarks) VALUES (?,?,?,?,?,?,?)", (cdr_type, cdrseq, pdbid_val, chain_val, cdr_len, tcr_type, remarks))
            
            
def update_db_last_updated():
    cur.execute('DROP TABLE IF EXISTS db_last_updated');
    cur.execute('''CREATE TABLE IF NOT EXISTS db_last_updated( last_updated_date text)''')
    cur.execute("INSERT INTO db_last_updated( last_updated_date ) VALUES (CURRENT_DATE)")
    return
    

conn = sqlite3.connect('TCR_database.db')
cur = conn.cursor()



print "updating tcr_chains table" 
#cur.execute('DROP TABLE IF EXISTS tcr_chains');
cur.execute('''CREATE TABLE IF NOT EXISTS tcr_chains( id integer primary key autoincrement, pdbid text, chain text, length integer, seq text, tcr_type text, desc text, release_date text, UNIQUE(pdbid, chain) )''')
update_tcr_chains(tcr_alpha_chain, "A", "Alpha" )
update_tcr_chains(tcr_beta_chain, "B", "Beta" )

print ""
print "updating cdr_data table" 
cur.execute('''CREATE TABLE IF NOT EXISTS cdr_data(id integer primary key autoincrement, cdr_type text, cdrseq text, pdbid text, chain text, cdr_length int, tcr_type text, remarks text)''')
update_cdr_data(tcr_alpha_chain, "A", "Alpha" )
update_cdr_data(tcr_beta_chain, "B", "Beta" )

print ""
print "updating germline_gene table" 
cur.execute('DROP TABLE IF EXISTS germline_gene');
cur.execute('''CREATE TABLE IF NOT EXISTS germline_gene(id integer primary key autoincrement, pdbid text, chain text, tcr_type text, subgroup_name text, gene_name text, allele_name text, species text, exact_match text)''')
#update_TRAV_imgt_gene_data()
update_germline_gene(tcr_alpha_chain, "A", trav_regex, "json_gene_files/TRAV_HUMAN.json", "json_gene_files/TRAV_MOUSE.json")
update_germline_gene(tcr_beta_chain, "B", trbv_regex, "json_gene_files/TRBV_HUMAN.json", "json_gene_files/TRBV_MOUSE.json")


update_db_last_updated()

print ""
print "updating complexes" 
update_complexes()

seq="GVTQTPRYLIKTRGQQVTLSCSPISGHRSVSWYQQTPGQGLQFLFEYFSETQRNKGNFPGRFSGRQFSNSRSEMNVSTLELGDSALYLCASSLEGGYYNEQFFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPA"
get_truncated_trv(seq);

conn.commit()
conn.close()

