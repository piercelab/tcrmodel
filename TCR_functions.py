import os
import gzip
import re
import string
from Bio.PDB import PDBParser, PDBIO, Select, CaPPBuilder, PPBuilder
from Bio import SeqIO
parser = PDBParser(PERMISSIVE=1,QUIET=True)
io = PDBIO()
ppb = CaPPBuilder()
from subprocess import Popen, PIPE
from Bio.PDB.Polypeptide import three_to_one
from shutil import copyfile
import scoring_matrices

PAM30 = ([6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2],
         [-7,8,-6,-10,-8,-2,-9,-9,-2,-5,-8,0,-4,-9,-4,-3,-6,-2,-10,-8],
         [-4,-6,8,2,-11,-3,-2,-3,0,-5,-7,-1,-9,-9,-6,0,-2,-8,-4,-8],
         [-3,-10,2,8,-14,-2,2,-3,-4,-7,-12,-4,-11,-15,-8,-4,-5,-15,-11,-8],
         [-6,-8,-11,-14,10,-14,-14,-9,-7,-6,-15,-14,-13,-13,-8,-3,-8,-15,-4,-6],
         [-4,-2,-3,-2,-14,8,1,-7,1,-8,-5,-3,-4,-13,-3,-5,-5,-13,-12,-7],
         [-2,-9,-2,2,-14,1,8,-4,-5,-5,-9,-4,-7,-14,-5,-4,-6,-17,-8,-6],
         [-2,-9,-3,-3,-9,-7,-4,6,-9,-11,-10,-7,-8,-9,-6,-2,-6,-15,-14,-5],
         [-7,-2,0,-4,-7,1,-5,-9,9,-9,-6,-6,-10,-6,-4,-6,-7,-7,-3,-6],
         [-5,-5,-5,-7,-6,-8,-5,-11,-9,8,-1,-6,-1,-2,-8,-7,-2,-14,-6,2],
         [-6,-8,-7,-12,-15,-5,-9,-10,-6,-1,7,-8,1,-3,-7,-8,-7,-6,-7,-2],
         [-7,0,-1,-4,-14,-3,-4,-7,-6,-6,-8,7,-2,-14,-6,-4,-3,-12,-9,-9],
         [-5,-4,-9,-11,-13,-4,-7,-8,-10,-1,1,-2,11,-4,-8,-5,-4,-13,-11,-1],
         [-8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8],
         [-2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-7,-6,-8,-10,8,-2,-4,-14,-13,-6],
         [0,-3,0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6],
         [-1,-6,-2,-5,-8,-5,-6,-6,-7,-2,-7,-3,-4,-9,-4,0,7,-13,-6,-3],
         [-13,-2,-8,-15,-15,-13,-17,-15,-7,-14,-6,-12,-13,-4,-14,-5,-13,13,-5,-15],
         [-8,-10,-4,-11,-4,-12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7],
         [-2,-8,-8,-8,-6,-7,-6,-5,-6,2,-2,-9,-1,-8,-6,-6,-3,-15,-7,7])


aa_map = {
    'A' : '0',
    'R' : '1',
    'N' : '2',
    'D' : '3',
    'C' : '4',
    'Q' : '5',
    'E' : '6',
    'G' : '7',
    'H' : '8',
    'I' : '9',
    'L' : '10',
    'K' : '11',
    'M' : '12',
    'F' : '13',
    'P' : '14',
    'S' : '15',
    'T' : '16',
    'W' : '17',
    'Y' : '18',
    'V' : '19'
    }


def seq_by_id_from_fasta_files_list(inp_id, fasta_files_list):
    for fasta_file in fasta_files_list:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            if rec.id == inp_id:
                return rec.seq

def pdbchain_to_fasta(inpdb, chainid):
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}
    pdbseq = ''
    prev = '-1'
    input_file = open(inpdb)
    for line in input_file:
        if line[:4] != 'ATOM': continue
        if line[21:22] != chainid: continue
        if line[22:27]!= prev:
            pdbseq += letters[line[17:20]]
        prev = line[22:27]
    return pdbseq

def pdb_to_fasta(inpdb):
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}
    pdbseq = ''
    prev = '-1'
    input_file = open(inpdb)
    for line in input_file:
        if line[:4] != 'ATOM': continue
        if line[22:27]!= prev:
            pdbseq += letters[line[17:20]]
        prev = line[22:27]
    return pdbseq

def make_pdb(pdb_path, chain_letters, outfile=None, overwrite=False, struct=None):
    chain_letters = [chain for chain in chain_letters]
    if outfile is None:
        outfile = "%s_%s.pdb" % (pdb_path, "".join(chain_letters))

    plural = "s" if (len(chain_letters) > 1) else ""  # for printing
    #print("Extracting chain%s %s from %s..." % (plural,", ".join(chain_letters), pdb_path))
    if struct is None:
        st = parser.get_structure("PDB", pdb_path)
        struct = st[0]
    io.set_structure(struct)
    #print "OUT PATH:",outfile
    io.save(outfile, select=SelectChains(chain_letters), write_end=False)
    
def get_pdb(pdbid,chainid=None,outfile=None,pdb_release_path=None):
    #pdb_release_path="/TCRmodeller/PDB_RELEASE/pdb_structures"
    gzpdbfile_path =  pdb_release_path + '/%s/pdb%s.ent.gz' %(pdbid[1:3], pdbid)
    gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
    st = parser.get_structure('PDB', gzpdbfile)
    struct = st[0]
    make_pdb(pdbid,chainid,outfile,True,struct)

class SelectChains(Select):
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters
    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)
    #to remove HETATM records
    def accept_residue(self,residue):
        if residue.id[0] == ' ':
            return 1
        else:
            return 0

def assign_CDRs_using_REGEX(seq, tag, nocap=False):
    import re
    tcra_regex_1 = "^[A-Z]{0,11}(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regex_2 = "^[A-Z]{0,11}(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_3 = "((^[A-Z]{1,30}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_4 = "((^[A-Z]{1,30}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_5 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]*))"
    tcra_regexList = [tcra_regex_1,tcra_regex_2]

    tcrb_regex_1 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C)([A-Z]C[A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    #for pdb 3to4
    tcrb_regex_2 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,31}DTELR)(L[A-Z]{12}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regex_3 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    #tcrb_regex_4 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{0,3}))[A-Z]*"
    #tcrb_regex_5 = "^[A-Z]{0,30}(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regexList = [tcrb_regex_1,tcrb_regex_2,tcrb_regex_3]

    if nocap is True:
        tcra_nocap_1 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]*))"
        tcra_nocap_2 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]*))"
        tcra_regexList = [tcra_nocap_1,tcra_nocap_2]

        tcrb_nocap_1 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C)([A-Z]C[A-Z]{1,32})(FG[A-Z]G[A-Z]*))"
        tcrb_nocap_2 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,31}DTELR)(L[A-Z]{12}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]*))"
        tcrb_nocap_3 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]*))"
        tcrb_regexList = [tcrb_nocap_1,tcrb_nocap_2, tcrb_nocap_3]

    if tag == 'A': tcr_regexList = tcra_regexList;
    if tag == 'B': tcr_regexList = tcrb_regexList;
    gotMatch = False
    for tcr_regex in tcr_regexList:
        #print tcr_regex 
        res = re.search(tcr_regex, str(seq))
        if res:
            #print "HIt"
            # we re.search again here with truncated string as input, this is to get span (location of match) from the truncated string 
            trunc = re.search(tcr_regex, str(res.group(1)))
            gotMatch = True
            return trunc

    if not gotMatch:        
        return None

def assign_CDRs_using_REGEX_webserver_version(seq, tag):
    import re
    tcra_regex_1 = "^[A-Z]*(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regex_2 = "^[A-Z]*(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regexList = [tcra_regex_1,tcra_regex_2]
    tcrb_regex_1 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C)([A-Z]C[A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    #for pdb 3to4
    tcrb_regex_2 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,31}DTELR)(L[A-Z]{12}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regex_3 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regexList = [tcrb_regex_1,tcrb_regex_2,tcrb_regex_3]
    if tag == 'A': tcr_regexList = tcra_regexList;
    if tag == 'B': tcr_regexList = tcrb_regexList;
    gotMatch = False
    for tcr_regex in tcr_regexList:
        res = re.search(tcr_regex, str(seq))
        if res:
            trunc = re.search(tcr_regex, str(res.group(1)))
            gotMatch = True
            return trunc
    if not gotMatch:        
        return None

def get_vdomain_using_regex(seq, tag):
    tcra_regex = "^[A-Z]*(([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y][A-Z][C])[A-Z]{1,32}([F|W][A-Z]{2}[G][A-Z]{6}))[A-Z]*"
    tcra_regexList = [tcra_regex]
    tcrb_regex = "^[A-Z]*(([A-Z]{19}[C][A-Z]{1,19}[W][A-Z]{1,70}[Y][A-Z][C])[A-Z]{1,32}([F|W][A-Z]{2}[G][A-Z]{6}))[A-Z]*"
    tcrb_regexList = [tcrb_regex]
    if tag == 'A': tcr_regexList = tcra_regexList;
    if tag == 'B': tcr_regexList = tcrb_regexList;
    gotMatch = False
    for tcr_regex in tcr_regexList:
        res = re.search(tcr_regex, str(seq))
        if res:
            trunc = re.search(tcr_regex, str(res.group(1)))
            gotMatch = True
            return trunc
    if not gotMatch:        
        return None
        

def assign_CDRs_using_REGEX_webserver_2version(seq, tag):
    import re
    tcra_regex_1 = "^[A-Z]*(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z])([C|W][A-Z]{1,32}[F|W])(G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regex_2 = "^[A-Z]*(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z])([C|W][A-Z]{1,32}[F|W])(G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regexList = [tcra_regex_1,tcra_regex_2]
    tcrb_regex_1 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12})(C[A-Z]C[A-Z]{1,32}F)(G[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    #for pdb 3to4
    tcrb_regex_2 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,31}DTELR)(L[A-Z]{12})(C[A-Z]{1,32}F)(G[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regex_3 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14})(C[A-Z]{1,32}F)(G[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regexList = [tcrb_regex_1,tcrb_regex_2,tcrb_regex_3]
    if tag == 'A': tcr_regexList = tcra_regexList;
    if tag == 'B': tcr_regexList = tcrb_regexList;
    gotMatch = False
    for tcr_regex in tcr_regexList:
        res = re.search(tcr_regex, str(seq))
        if res:
            trunc = re.search(tcr_regex, str(res.group(1)))
            gotMatch = True
            return trunc
    if not gotMatch:        
        return None


def run_anarci(query, anarci_program=None, outfile=None, restrict=None):
    if anarci_program is None:
        anarci_program = "ANARCI"
    if outfile is None:
        outfile = "anarci.out"
    if os.path.isfile(outfile):
        os.remove(outfile)
    if restrict is None:
        process = Popen([anarci_program, "-i" , str(query), "-o" , str(outfile) ,  "-s" , "a" ], stdout=PIPE, stderr=PIPE)
    else:
        process = Popen([anarci_program, "-i" , str(query), "-o" , str(outfile) ,  "-s" , "a" , "-r", str(restrict)], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

def get_seq_from_pdb(inpdb):
    ppb = PPBuilder()# Using C-N
    #ppb = CaPPBuilder()# Using CA-CA
    inpstruct = parser.get_structure('TCR', inpdb)
    seqlist = []
    for chaininfo in inpstruct.get_chains():
        s = inpstruct[0][chaininfo.id]
        seq = ""
        for ppe in ppb.build_peptides(s):
            seq += str(ppe.get_sequence())
        seqlist.append(seq)
    return seqlist


def get_tcr_vdomain(record, TAG):
    pdbfile = record.id +".pdb"
    anarci_out_file = "anarci.out"
    get_pdb(record.id[:4], record.id[5:6], pdbfile, "/TCRmodeller/PDB_RELEASE/pdb_structures" )
    #inseq = get_seq_from_pdb(pdbfile)
    #pdbseq = inseq[0] 
    pdbseq = pdb_to_fasta(pdbfile)
    run_anarci(pdbseq,None,anarci_out_file,None)
    ahooutfile = record.id + "." + TAG + ".aho.pdb"
    parse_anarci_output(anarci_out_file, pdbfile, pdbseq, TAG, ahooutfile)


def calc_alphacys_betacys_distance(apdb,bpdb):
    achainid = apdb[5:6]
    bchainid = bpdb[5:6]
    astructure = parser.get_structure("Achain", apdb) 
    bstructure = parser.get_structure("Bchain", bpdb) 
    aatom = astructure[0][achainid][106]['CA']
    batom = bstructure[0][bchainid][106]['CA']
    return aatom - batom


def make_complex(apdb,bpdb,outfile=None):
    achainid = apdb[5:6]
    bchainid = bpdb[5:6]
    astructure = parser.get_structure("Achain", apdb) 
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
    if outfile is None:
        #outfile = apdb[:4]+"_"+achainid+bchainid+"_aho.pdb"
        outfile = apdb[:4]+"_"+achainid+bchainid+"_aho.pdb"
    io.save(outfile)
    return outfile

def cap_pdb(pdbfile,capbegin,capend,outfile=None):
    st = parser.get_structure('PDB', pdbfile)
    i = 0#pdb_counter
    for chain in st[0]:
        for residue in list(chain):
            if (i < int(capbegin)) or (i >= int(capend)):
                chain.detach_child(residue.id)
            i += 1

    io.set_structure(st)
    if outfile is None:
        outfile = pdbfile+"_cap.pdb"
    io.save(outfile)


def renumber_pdbfile_to_aho(pdbfile, TAG=None, ahooutfile=None, extend=False):
    pdbseq = pdb_to_fasta(pdbfile)
    filetag = os.path.basename(pdbfile)
    anarci_out_file = filetag+"_anarci.out"
    run_anarci(pdbseq,"ANARCI",anarci_out_file,None)
    aho_begin = ''
    aho_end = ''
    dom = "~!@#$%^&*()_+"#just random initialization 
    anarcilist = []
    FoundDomain = None
    if TAG is None:
        TAG = "[A-Z]"
    with open(anarci_out_file) as f:
        regex_for_header = "\#\|.*\|("+TAG+")\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
        for line in f:
            res = re.search(regex_for_header, line)
            if res:
                dom = res.groups()[0]
                aho_begin = res.groups()[1]
                aho_end = res.groups()[2]
                FoundDomain = True
                break
        if not FoundDomain:
            print "Failed to identify TCR domain with ANARCI program: ", pdbfile
            return None
            
    with open(anarci_out_file) as f:
        for line in f:
            if line.startswith("#"):continue
            if line.startswith("//"):continue
            if line.startswith(dom):
                if line[10] == "-":continue
                values = [line[0],line[2:7],line[8],line[10]]
                anarcilist.append(values)
    #print tcrdomain
    #print anarcilist
    #print anarcilist[0][1], anarcilist[-1][1]

    st = parser.get_structure('PDB', pdbfile)
    
    if extend is True:
        counter = int(aho_begin)
        while (counter > 0):
            counter -= 1
            aacode = pdbseq[counter]
            aanum = int(anarcilist[0][1].strip()) - 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.insert(0, values)
        counter = int(aho_end)
        while (counter < (len(pdbseq)-1)):
            counter += 1
            aacode = pdbseq[counter]
            aanum = int(anarcilist[-1][1].strip()) + 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.append(values)
        j = 0#anarci_counter
        for chain in st[0]:
            for residue in list(chain):
                #print i, j, three_to_one(residue.resname), anarcilist[j], residue, residue.resname, three_to_one(residue.resname)
                assert (three_to_one(residue.resname) == anarcilist[j][3]), "Residues does not match %r %r %r" %(j, three_to_one(residue.resname), anarcilist[j][3])
                residue.id = (' ',int(anarcilist[j][1].strip()), anarcilist[j][2])
                j += 1

    if extend is False:
        i = 0#pdb_counter
        j = 0#anarci_counter
        for chain in st[0]:
            for residue in list(chain):
                if (i < int(aho_begin)) or (i > int(aho_end)):
                    chain.detach_child(residue.id)
                    i += 1
                    continue
                else:
                    #print i, j, three_to_one(residue.resname), anarcilist[j], residue, residue.resname, three_to_one(residue.resname)
                    assert (three_to_one(residue.resname) == anarcilist[j][3]), "Residues does not match %r %r %r" %(j, three_to_one(residue.resname), anarcilist[j][3])
                    residue.id = (' ',int(anarcilist[j][1].strip()), anarcilist[j][2])
                    i += 1
                    j += 1


    io.set_structure(st)
    if ahooutfile is None:
        ahooutfile = filetag+".aho.pdb"
    io.save(ahooutfile, write_end=False)
    return ahooutfile

def check_tcr_alpha_beta_domain(inseq):
    anarci_out_file = "anarci.out"
    run_anarci(inseq,ANARCI,anarci_out_file,None)
    a_regex = "\#\|.*\|(A)\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    b_regex = "\#\|.*\|(B)\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    Found_A = False
    Found_B = False
    with open(anarci_out_file) as f:
        for line in f:
            a_res = re.search(a_regex, line)
            b_res = re.search(b_regex, line)
            if a_res:
                Found_A = True
            if b_res:
                Found_B = True
    if ( (Found_A is True) or (Found_B is True) ):
        return True
    else:
        return False


def get_tcr_domain_seq(inseq, TAG):
    anarci_out_file = "anarci.out"
    run_anarci(inseq,ANARCI,anarci_out_file,None)
    regex_for_header = "\#\|.*\|("+TAG+")\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
    with open(anarci_out_file) as f:
        for line in f:
            res = re.search(regex_for_header, line)
            if res:
                dom = res.groups()[0]
                aho_begin = int(res.groups()[1])
                aho_end = int(res.groups()[2])
                break

    tcrdomain=inseq[aho_begin:aho_end+1]
    return tcrdomain

def change_chain_id(pdbfile,old_chainid,new_chainid):
    IN = open(pdbfile)
    outfilename = "new_"+pdbfile
    OUT = open(outfilename, 'w+')

    for line in IN.readlines():
        if line[0:4] == 'ATOM':
            if line[21:22] == old_chainid:
                line  = line[:21] + new_chainid + line[22:]
        OUT.write(line)
    return outfilename



def assign_CDRs_using_REGEX_shorter_aho(seq, tag, nocap=False):
    import re
    tcra_regex_1 = "^[A-Z]{0,11}(([A-Z]{19}C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W][A-Z]{2})([A-Z]{1,29})([A-Z][F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    tcra_regex_2 = "^[A-Z]{0,11}(([A-Z]{19}C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W][A-Z]{2})([A-Z]{1,29})([A-Z][F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_3 = "((^[A-Z]{1,30}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_4 = "((^[A-Z]{1,30}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*"
    #tcra_regex_5 = "((^[A-Z]*C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]*))"
    tcra_regexList = [tcra_regex_1,tcra_regex_2]

    tcrb_regex_1 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C[A-Z]C)([A-Z]{1,29})([A-Z]FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regex_2 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C[A-Z]{2})([A-Z]{1,29})([A-Z]FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    #tcrb_regex_3 = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{0,3}))[A-Z]*"
    #tcrb_regex_4 = "^[A-Z]{0,30}(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*"
    tcrb_regexList = [tcrb_regex_1,tcrb_regex_2]

    if nocap is True:
        tcra_nocap_1 = "((^[A-Z]*C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W][A-Z]{2})([A-Z]{1,29})([A-Z][F|W]G[A-Z]G[A-Z]*))"
        tcra_nocap_2 = "((^[A-Z]*C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W][A-Z]{2})([A-Z]{1,29})([A-Z][F|W]G[A-Z]G[A-Z]*))"
        tcra_regexList = [tcra_nocap_1,tcra_nocap_2]

        tcrb_nocap_1 = "((^[A-Z]*C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C[A-Z]C)([A-Z]{1,29})([A-Z]FG[A-Z]G[A-Z]*))"
        tcrb_nocap_2 = "((^[A-Z]*C[A-Z])([A-Z]{1,16})([A-Z]{2}W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C[A-Z]{2})([A-Z]{1,29})([A-Z]FG[A-Z]G[A-Z]*))"
        tcrb_regexList = [tcrb_nocap_1,tcrb_nocap_2]

    if tag == 'A': tcr_regexList = tcra_regexList;
    if tag == 'B': tcr_regexList = tcrb_regexList;
    gotMatch = False
    for tcr_regex in tcr_regexList:
        #print tcr_regex 
        #print seq

        res = re.search(tcr_regex, str(seq))
        if res:
            #print "HIt"
            # we re.search again here with truncated string as input, this is to get span (location of match) from the truncated string 
            trunc = re.search(tcr_regex, str(res.group(1)))
            gotMatch = True
            return trunc

    if not gotMatch:        
        return None

def score_alignment_from_fasta_files_list(inseq, fasta_files_list, inp_sc_mat=None, pdb_blacklist=None):
    best_score = -99999
    best_record = ""
    for fasta_file in fasta_files_list:
        curr_score, curr_record = score_alignment_from_fasta_file(inseq, fasta_file, inp_sc_mat=None, pdb_blacklist=pdb_blacklist)
        if (curr_score > best_score):
            best_score = curr_score
            best_record = curr_record
    return best_score, best_record


def score_alignment_from_fasta_file(inseq, fasta_file, inp_sc_mat=None, pdb_blacklist=None):
    best_score = -99999
    best_record = ""
    score = -99999

    handle = open(fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        if ( len(inseq) != len(record.seq) ): continue
        if pdb_blacklist:
            pdb_ignorelist = [x.lower() for x in pdb_blacklist]
            if (record.id[:4].lower() in pdb_ignorelist): continue
        #score = calculate_identity_score(inseq, record.seq)
        score = score_alignment(inseq, record.seq)
        if (score > best_score):
            best_score = score
            best_record = record
    return best_score, best_record

def score_alignment(inseq, dbseq, inp_sc_mat=None):
    assert (len(inseq) == len(dbseq)), "sequence length not same: %s, %s" % (inseq, dbseq)
    if inp_sc_mat is None:
        pssm = scoring_matrices.PAM30
    else:
        pssm = getattr(scoring_matrices,inp_sc_mat)
    score = 0
    for x in xrange(0, len(inseq)):
        score += pssm[int(aa_map[inseq[x:x+1]])][int(aa_map[dbseq[x:x+1]])]
    return score

def seq_match_from_fasta_files_list(inseq, fasta_files_list, pdb_blacklist=None):
    for fasta_file in fasta_files_list:
        curr_record = seq_match_from_fasta_file(inseq, fasta_file, pdb_blacklist=pdb_blacklist)
        if curr_record is not None:
            return curr_record
    return curr_record

def seq_match_from_fasta_file(inseq, fasta_file, pdb_blacklist=None):
    best_record = None
    handle = open(fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        if ( len(inseq) != len(record.seq) ): continue
        if pdb_blacklist:
            pdb_ignorelist = [x.lower() for x in pdb_blacklist]
            if (record.id[:4].lower() in pdb_ignorelist): continue
        if (str(inseq) == str(record.seq)):
            best_record = record
            break
    return best_record


def calculate_identity_score( query, content):
    if (len(query) != len(content)): return -99999
    percent_identity = 0
    num_simi_res = 0
    for x in xrange(0, len(query)):
        if ( query[x] == content[x] ): num_simi_res += 1
    percent_identity = ( num_simi_res / float(len(query)) ) * 100;
    return percent_identity;

def get_cdr_from_seq_by_aho_num(inseq, TAG=None, extend=None):
    if TAG == "A":
        cdr_aho_num = [24,42,56,71,107,138,81,90,56,90]
    elif TAG == "B":
        cdr_aho_num = [24,42,56,70,107,138,81,90,56,90]
    inseq = str(inseq.rstrip())
    anarci_out_file = "cdranarci.out"
    run_anarci(inseq, None, anarci_out_file,None)
    ahonumlist = []
    count = 0
    #alpha                                                                                                          
    cdr1_begin_pos = None
    cdr1_end_pos = None
    cdr2_begin_pos = None
    cdr2_end_pos = None
    cdr3_begin_pos = None
    cdr3_end_pos = None
    hv4_begin_pos = None
    hv4_end_pos = None
    cdr2hv4_begin_pos = None
    cdr2hv4_end_pos = None
    dom = ''
    aho_begin = ''
    aho_end = ''
    anytag = "[A-Z]"
    #if TAG is None:
        #TAG = "[A-Z]"
    with open(anarci_out_file) as f:
        regex_for_header = "\#\|.*\|("+anytag+")\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
        for line in f:
            res = re.search(regex_for_header, line)
            if res:
                dom = res.groups()[0]
                aho_begin = res.groups()[1]
                aho_end = res.groups()[2]
                break

    anarcilist = []
    count = 0
    with open(anarci_out_file) as f:
        for line in f:
            if line.startswith("#"):continue
            if line.startswith("//"):continue
            if line.startswith(dom):
                if line[10] == "-":continue
                values = [line[0],line[2:7],line[8],line[10]]
                if int(values[1].strip()) == cdr_aho_num[0]: cdr1_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[1]: cdr1_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[2]: cdr2_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[3]: cdr2_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[8]: cdr2hv4_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[9]: cdr2hv4_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[6]: hv4_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[7]: hv4_end_pos = count
                if int(values[1].strip()) == cdr_aho_num[4]: cdr3_begin_pos = count
                if int(values[1].strip()) == cdr_aho_num[5]: cdr3_end_pos = count
                count += 1
                anarcilist.append(values)
            else:
                return None
    if extend is True:
        cdr1_begin_pos += int(aho_begin)
        cdr1_end_pos += int(aho_begin)
        cdr2_begin_pos += int(aho_begin)
        cdr2_end_pos += int(aho_begin)
        cdr3_begin_pos += int(aho_begin)
        cdr3_end_pos += int(aho_begin)
        hv4_begin_pos += int(aho_begin)
        hv4_end_pos += int(aho_begin)
        cdr2hv4_begin_pos += int(aho_begin)
        cdr2hv4_end_pos += int(aho_begin)
    
        counter = int(aho_begin)
        while (counter > 0):
            counter -= 1
            aacode = inseq[counter]
            aanum = int(anarcilist[0][1].strip()) - 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.insert(0, values)
        counter = int(aho_end)
        while (counter < (len(inseq)-1)):
            counter += 1
            aacode = inseq[counter]
            aanum = int(anarcilist[-1][1].strip()) + 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.append(values)

    cdr1 = ''
    cdr2 = ''
    cdr2hv4 = ''
    hv4 = ''
    cdr3 = ''
    fw = ''
    cdr3_extnd = ''
    fw1 = ''
    fw2 = ''
    fw3 = ''
    fw4 = ''

    for item in anarcilist[:cdr1_begin_pos:]:
        fw1 += item[3].rstrip()
        fw += item[3].rstrip()
    for item in anarcilist[cdr1_end_pos+1:cdr2hv4_begin_pos]:
        fw2 += item[3].rstrip()
        fw += item[3].rstrip()
    for item in anarcilist[cdr2hv4_end_pos+1:cdr3_begin_pos]:
        fw3 += item[3].rstrip()
        fw += item[3].rstrip()
    for item in anarcilist[cdr3_end_pos+1:]:
        fw4 += item[3].rstrip()
        fw += item[3].rstrip()

    for item in anarcilist[cdr1_begin_pos:cdr1_end_pos+1]:
        cdr1 += item[3].rstrip()
    for item in anarcilist[cdr2_begin_pos:cdr2_end_pos+1]:
        cdr2 += item[3].rstrip()
    for item in anarcilist[cdr2hv4_begin_pos:cdr2hv4_end_pos+1]:
        cdr2hv4 += item[3].rstrip()
    for item in anarcilist[hv4_begin_pos:hv4_end_pos+1]:
        hv4 += item[3].rstrip()
    for item in anarcilist[cdr3_begin_pos:cdr3_end_pos+1]:
        cdr3 += item[3].rstrip()
    for item in anarcilist[cdr3_begin_pos-1:cdr3_end_pos+1+1]:
        cdr3_extnd += item[3].rstrip()
    return cdr1,cdr2,cdr2hv4,hv4,cdr3,fw,cdr3_extnd,fw1,fw2,fw3,fw4


def extract_deposition_date_from_pdb(record_id):
    pdbid = record_id[:4]
    if (pdbid == "4zez"):
        pdbid = "5jzi"
    chainid = pdbid[5:6]
    pdb_release_path = "/TCRmodeller/PDB_RELEASE/pdb_structures"
    gzpdbfile_path =  pdb_release_path + '/%s/pdb%s.ent.gz' %(pdbid[1:3], pdbid)
    if os.path.isfile(gzpdbfile_path):
        gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
        st = parser.get_structure(pdbid, gzpdbfile)
        deposition_date = st.header['deposition_date']
        resolution = st.header['resolution']
        return deposition_date, resolution
    else:
        return None

def extract_resolution_from_pdb(record_id):
    pdbid = record_id[:4]
    if (pdbid == "4zez"):
        pdbid = "5jzi"
    chainid = pdbid[5:6]
    pdb_release_path = "/TCRmodeller/PDB_RELEASE/pdb_structures"
    gzpdbfile_path =  pdb_release_path + '/%s/pdb%s.ent.gz' %(pdbid[1:3], pdbid)
    if os.path.isfile(gzpdbfile_path):
        gzpdbfile = gzip.open(gzpdbfile_path, 'rb')
        st = parser.get_structure(pdbid, gzpdbfile)
        resolution = st.header['resolution']
        return resolution
    else:
        return None


def find_orientation_template(seqa, seqb, orientation_template_file, inp_sc_mat=None, pdb_blacklist=None):
    identity_cutoff = 100
    best_score = -9999999
    best_alpha_tag = ''
    best_beta_tag = ''
    #print seqa, seqb
    with open(orientation_template_file, 'r') as f:
        for line in f:
            if not line.strip(): continue
            linelist = line.split()
            #if ( (seqa == linelist[2]) and (seqb == linelist[3]) ): continue
            if ( (len(seqa) == len(linelist[2])) and (len(seqb) == len(linelist[3])) ):
                if pdb_blacklist:
                    pdb_ignorelist = [x.lower() for x in pdb_blacklist]
                    if (linelist[0][:4].lower() in pdb_ignorelist): continue
                    if (linelist[1][:4].lower() in pdb_ignorelist): continue
                identity_score = ( calculate_identity_score( str(seqa+seqb), str(linelist[2]+linelist[3]) ) )
                if (calculate_identity_score( str(seqa+seqb), str(linelist[2]+linelist[3]) ) > identity_cutoff ): continue
                score = score_alignment(seqa,linelist[2],inp_sc_mat) + score_alignment(seqb,linelist[3],inp_sc_mat)
                if (score > best_score):
                    best_score = score
                    best_alpha_tag = linelist[0]
                    best_beta_tag = linelist[1]
    return best_alpha_tag,best_beta_tag


def rename_pdb_chains(cur_chains, new_chains, infile, outfile=None):
    assert len(cur_chains) == len(new_chains), "No, of current/new chains do not match."
    tmpoutfile = "chain_rename_tmp.pdb"
    IN = open(infile)
    OUT = open(tmpoutfile, 'w+')
    for line in IN.readlines():
        if line[0:4] == 'ATOM' or line[0:4] == 'HETATM':            
            newline = line
            for count, value in enumerate(cur_chains):
                if line[21:22] == cur_chains[count]:
                    newline  = line[:21] + new_chains[count] + line[22:]
                    break
            OUT.write(newline)
        else:
            OUT.write(line)
    IN.close()
    OUT.close()
    if outfile is None:
        copyfile(tmpoutfile, infile)
    else:
        copyfile(tmpoutfile, outfile)
    os.remove(tmpoutfile)
    return

def renumber_pdbchain_to_aho(pdbfile, chainid=None, TAG=None, ahooutfile=None, extend=False):
    if chainid is None:
        pdbseq = pdb_to_fasta(pdbfile)
    else:
        pdbseq = pdbchain_to_fasta(pdbfile,chainid)

    filetag = os.path.basename(pdbfile)
    anarci_out_file = filetag+"_anarci.out"
    run_anarci(pdbseq,"ANARCI",anarci_out_file,None)
    aho_begin = ''
    aho_end = ''
    dom = "~!@#$%^&*()_+"#just random initialization 
    anarcilist = []
    FoundDomain = None
    if TAG is None:
        TAG = "[A-Z]"
    with open(anarci_out_file) as f:
        regex_for_header = "\#\|.*\|("+TAG+")\|.*\|.*\|([0-9]+)\|([0-9]+)\|"
        for line in f:
            res = re.search(regex_for_header, line)
            if res:
                dom = res.groups()[0]
                aho_begin = res.groups()[1]
                aho_end = res.groups()[2]
                FoundDomain = True
                break
        if not FoundDomain:
            print "Failed to identify TCR domain using ANARCI program: ", pdbfile
            return None
    #Kappa and lambda are considered together in ANARCI
    if dom == "K":
        dom = "L"

    with open(anarci_out_file) as f:
        for line in f:
            if line.startswith("#"):continue
            if line.startswith("//"):continue
            if line.startswith(dom):
                if line[10] == "-":continue
                values = [line[0],line[2:7],line[8],line[10]]
                anarcilist.append(values)
    #print dom
    #print anarcilist
    #print anarcilist[0][1], anarcilist[-1][1]

    st = parser.get_structure('PDB', pdbfile)
    if chainid is None:
        chain = st[0][0]#first chain
    else:
        chain = st[0][chainid]

    if extend is True:
        counter = int(aho_begin)
        while (counter > 0):
            counter -= 1
            aacode = pdbseq[counter]
            aanum = int(anarcilist[0][1].strip()) - 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.insert(0, values)
        counter = int(aho_end)
        while (counter < (len(pdbseq)-1)):
            counter += 1
            aacode = pdbseq[counter]
            aanum = int(anarcilist[-1][1].strip()) + 1
            values = [TAG,str(aanum)," ",aacode]      
            anarcilist.append(values)
        j = 0#anarci_counter
        #for chain in st[0]:
        for residue in list(chain):
            #print i, j, three_to_one(residue.resname), anarcilist[j], residue, residue.resname, three_to_one(residue.resname)
            assert (three_to_one(residue.resname) == anarcilist[j][3]), "Residues does not match %r %r %r" %(j, three_to_one(residue.resname), anarcilist[j][3])
            residue.id = (' ',int(anarcilist[j][1].strip()), anarcilist[j][2])
            j += 1

    if extend is False:
        i = 0#pdb_counter
        j = 0#anarci_counter
        #for chain in st[0]:
        for residue in list(chain):
            if (i < int(aho_begin)) or (i > int(aho_end)):
                chain.detach_child(residue.id)
                i += 1
                continue
            else:
                print i, j, three_to_one(residue.resname), anarcilist[j], residue, residue.resname, three_to_one(residue.resname)
                assert (three_to_one(residue.resname) == anarcilist[j][3]), "Residues does not match %r %r %r" %(j, three_to_one(residue.resname), anarcilist[j][3])
                residue.id = (' ',int(anarcilist[j][1].strip()), anarcilist[j][2])
                i += 1
                j += 1

    io.set_structure(st)
    if ahooutfile is None:
        ahooutfile = filetag+".aho.pdb"
    io.save(ahooutfile)
    return ahooutfile
