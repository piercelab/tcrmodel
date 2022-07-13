from flask import Flask, g, render_template, request, redirect, url_for
from flask import json
import os
import paramiko
import string, datetime, random
import subprocess
import commands
import json
import glob
import shutil
import zipfile
from werkzeug.utils import secure_filename
from TCR_functions import *
import sys
import pwd
import grp
import re

app = Flask(__name__)

tcrdata = {}


rundir = "tcrmodel_runs"
#Works if rundirpath is one directory down the app (tcrmodel) folder
rundir_path = os.path.abspath(os.path.join(app.root_path, '..', rundir))
#rundir server path
rundir_spath = "/"+rundir

#template db path for pdb and seq
tmplt_db_pdb = "/www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr/pdb/"
tmplt_db_seq = "/www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr/seq/"

tcrdata['rundir_spath'] = rundir_spath


def get_unique_name(basename):
   random_tag = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(3))
   datetime_tag = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
   uniquedirname = "_".join([("".join([basename, random_tag])), datetime_tag])
   return uniquedirname

def create_and_cd_to_unique_dir(basename):
   uniquedirname = get_unique_name(basename)
   rundir = os.path.join(rundir_path,uniquedirname)
   if not os.path.exists(rundir):
      os.makedirs(rundir)
   os.chmod(rundir, 0777)
   os.chdir(rundir)
   return uniquedirname

def get_core_binding_pepseq(pseq):
   #subprocess.call(["/www/cgi-bin/netMHCIIpan-4.1/netMHCIIpan -inptype 1 -f ", pepfile], stdout=outfile)
   pepfile = "inpepseq.txt"
   with open(pepfile, 'w') as pfile:
      pfile.write(pseq)
   commandline = "/www/cgi-bin/netMHCIIpan-4.1/netMHCIIpan -inptype 1 -f " + pepfile 
   p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out,err = p.communicate()
   for line in out.splitlines():
      #skip blank lines, line start with #, - , Pos, Number
      if len(line.strip()) == 0 :
         continue
      if line.startswith("#"):
         continue
      if line.startswith("-"):
         continue
      if line.startswith(" Pos"):
         continue
      if line.startswith("Number"):
         continue
      core_pseq = line.split()[4]
   return core_pseq

def find_templates(inseq, multidb, tag1, tag2, inp_sc_mat=None, pdb_blacklist=None):         
   tmplt = score_alignment_from_fasta_files_list( inseq, multidb, None, pdb_blacklist )
   if all(tmplt):
      tcrdata[tag1+'_tmplt_id'] = str(tmplt[1].id)
      tcrdata[tag1+'_tmplt_pdb'] = str(tmplt[1].id)[:6].upper()
      tcrdata[tag1+'_tmplt_pdbid'] = str(tmplt[1].id)[:4].upper()
      tcrdata[tag1+'_tmplt_pdb_chain'] = str(tmplt[1].id)[5:6].upper()
      tcrdata[tag1+'_tmplt_seq'] = str(tmplt[1].seq)
   else:
      errormsg = "No Template found for TCR "+tag2+" segment"
      return render_template("error.html", errormsg=errormsg)

def submit_job_on_cluster(uniquejobid,jobscriptfile,ig_file):
   HOST="cluster.ibbr.umd.edu"
   USER="tcrmodel"
   PASSWD="o8!St1Jfqx5xo=x^8FXv%$_sn4W#BO"
   remote_dir = "/home/gowthamanr/TCRmodeller/"
   remote_location = os.path.join(remote_dir, uniquejobid)
   ssh_client = paramiko.SSHClient()
   ssh_client.load_system_host_keys()
   ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
   ssh_client.connect(HOST, username=USER, password=PASSWD)
   sftp = ssh_client.open_sftp()
   parent = os.path.join(rundir_path, uniquejobid)
   sftp.mkdir(remote_location)
   sftp.chdir(remote_location)
   sftp.put(jobscriptfile, jobscriptfile)
   sftp.put(os.path.join(parent, jobscriptfile), os.path.join(remote_location, jobscriptfile))
   ignore_file = os.path.join(parent, ig_file)
   if os.path.isfile( ignore_file ):
      sftp.put(ignore_file, os.path.join(remote_location, ig_file))
   command = "cd %s && sbatch %s" % (remote_location, jobscriptfile)
   ssh_stdin, ssh_stdout, ssh_stderr = ssh_client.exec_command(command)
   return

def submit_batchjob_on_cluster(uniquejobid,jobscriptfile,ig_file):
   HOST="cluster.ibbr.umd.edu"
   USER="tcrmodel"
   PASSWD="o8!St1Jfqx5xo=x^8FXv%$_sn4W#BO"
   remote_dir = "/home/gowthamanr/TCRmodeller/"
   remote_location = os.path.join(remote_dir, uniquejobid)
   ssh_client = paramiko.SSHClient()
   ssh_client.load_system_host_keys()
   ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
   ssh_client.connect(HOST, username=USER, password=PASSWD)
   sftp = ssh_client.open_sftp()
   parent = os.path.join(rundir_path, uniquejobid)
   for dirpath, dirnames, filenames in os.walk(parent):
      remote_path = os.path.join(remote_location, dirpath[len(parent)+1:])
      sftp.mkdir(remote_path)
      for filename in filenames:
         sftp.put(os.path.join(dirpath, filename), os.path.join(remote_path, filename))
   command = "cd %s && csh %s" % (remote_location, jobscriptfile)
   ssh_stdin, ssh_stdout, ssh_stderr = ssh_client.exec_command(command)
   return

def get_seq_from_genename(gene, trfile):

   #(cleanup genenames: Ex. TRBV06-05 to TRBV6-5)
   gene = re.sub('-0','-',gene)
   gene = re.sub('TRAV0','TRAV',gene)
   gene = re.sub('TRAJ0','TRAJ',gene)
   gene = re.sub('TRBV0','TRBV',gene)
   gene = re.sub('TRBJ0','TRBJ',gene)
   trseq = ''
   if trfile:
      with open(trfile, 'r') as f:
         trfiledata = json.loads(f.read());
         Found = False
         for item in trfiledata:
            if item.get('allele_name') == gene:
               trseq = item.get('seq')
               Found = True
               break;
         if not Found:
            for item in trfiledata:
               if item.get('gene_name') == gene:
                  trseq = item.get('seq')
                  Found = True
                  break;
         if not Found:
            #cleanup gene name
            if '*' not in gene:
               gene += "*01"
            gene = re.sub('-.+\*','*',gene)
            for item in trfiledata:
               if item.get('allele_name') == gene:
                  trseq = item.get('seq')
                  Found = True
                  break;
         if not Found:
            gene = re.sub('\*.+','',gene)
            for item in trfiledata:
               if item.get('gene_name') == gene:
                  trseq = item.get('seq')
                  Found = True
                  break;
         if not Found:
            #cleanup gene name(TRAV14-1*01 to TRAV14)
            gene = re.sub('-.+','',gene)
            for item in trfiledata:
               if item.get('subgroup_name') == gene:
                  trseq = item.get('seq')
                  Found = True
                  break;
   return trseq

def find_geneseq_from_cdr3seq(cdr3seq, trfile):
   trseq = ''
   if trfile:
      with open(trfile, 'r') as f:
         trfiledata = json.loads(f.read());
         for i in range(0, len(cdr3seq)):
            inseq =  cdr3seq[i:]
            Found = False
            for item in trfiledata:
               if inseq in item.get('halfseq'):
                  trseq = item.get('seq')
                  Found = True
                  break;
            if Found:
               break
   return trseq
   
@app.route('/external_submit/<path:genelist>')
def external_submit(genelist):
   glist = genelist.split(",")
   trav = glist[0]
   cdr3a = glist[1]
   traj = glist[2]
   trbv = glist[3]
   cdr3b = glist[4]
   trbj = glist[5]

   travfile = app.root_path + '/static/genemapdata/new_TRAV_HUMAN.json';
   travseq = get_seq_from_genename(trav, travfile)

   trajfile = app.root_path + '/static/genemapdata/new_TRAJ_HUMAN.json';
   if traj:
      trajseq = get_seq_from_genename(traj, trajfile)      
   if not traj or trajseq == '':
      trajseq = find_geneseq_from_cdr3seq(cdr3a, trajfile)


   trbvfile = app.root_path + '/static/genemapdata/new_TRBV_HUMAN.json';
   trbvseq = get_seq_from_genename(trbv, trbvfile)

   trbjfile = app.root_path + '/static/genemapdata/new_TRBJ_HUMAN.json';
   if trbj:
      trbjseq = get_seq_from_genename(trbj, trbjfile)      
   if not trbj or trbjseq == '':
      trbjseq = find_geneseq_from_cdr3seq(cdr3b, trbjfile)

   aseq = travseq+cdr3a+trajseq
   bseq = trbvseq+cdr3b+trbjseq

   return render_template("external_submit.html",aseq=aseq,bseq=bseq)

@app.route('/external_submit1/<trav>/<traj>/<cdr3a>/<trbv>/<trbj>/<cdr3b>')
def external_submit1(trav,traj,cdr3a,trbv,trbj,cdr3b):
   return render_template("external_submit.html",trav=trav,traj=traj,cdr3a=cdr3a,trbv=trbv,trbj=trbj,cdr3b=cdr3b)

'''
@app.route('/test')
def test():
   return render_template("test.html")
'''

@app.route('/')
@app.route('/index')
def index():
   return render_template("index.html")
   
@app.route('/about')
def about():
   return render_template("about.html")

@app.route('/help')
def help():
   return render_template("help.html")

@app.route('/links')
def links():
   return render_template("links.html")

@app.route('/processjob/<aseq>/<bseq>/<loopref_checked>/<pdb_blacklist>')
def processjob(aseq,bseq,loopref_checked,pdb_blacklist):
   aseq=''.join(i for i in aseq if i.isalpha())
   bseq=''.join(i for i in bseq if i.isalpha())
   tcrdata['aseq_user']=aseq
   tcrdata['bseq_user']=bseq
 
   a_regexp_res = get_vdomain_using_regex(aseq,"A") 
   if (not a_regexp_res):
      errormsg = "No TCR Variable domain sequence identified from the Alpha chain input! If you used 'Select from gene', the CDR3 sequence was invalid."
      return render_template("error.html", errormsg=errormsg)
   a_vd_seq = a_regexp_res.groups()[0]
   a_gm_seq = a_regexp_res.groups()[1]+a_regexp_res.groups()[2]

   b_regexp_res = get_vdomain_using_regex(bseq,"B") 
   if (not b_regexp_res):
      errormsg = "No TCR Variable domain sequence identified from the Beta chain input! If you used 'Select from gene', the CDR3 sequence was invalid."
      return render_template("error.html", errormsg=errormsg)
   b_vd_seq = b_regexp_res.groups()[0]
   b_gm_seq = b_regexp_res.groups()[1]+b_regexp_res.groups()[2]

   uniquedirname = create_and_cd_to_unique_dir("TCRS")#TCRS for single process
   tcrdata['jobid'] = uniquedirname

   #rtcrcommand = "-mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s -template_similarity_cutoff %s "  % (aseq,bseq,simil_cutoff) 
   rtcrcommand = "-mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s -minimize_model "  % (aseq,bseq) 

   if pdb_blacklist:
      pdb_blacklist = pdb_blacklist.split(',')
      with open("ignore_list.txt", 'w') as ignore_list_file:
         for ignorepdb in pdb_blacklist:
            ignore_list_file.write(ignorepdb+'\n')
            rtcrcommand += "-ignore_list ignore_list.txt "

   looprem_checked = "no";
   if (looprem_checked == "yes"):
      rtcrcommand += "-remodel_tcr_cdr3_loops "
   if (loopref_checked == "yes"):
      rtcrcommand += "-refine_tcr_cdr3_loops -loops::max_inner_cycles 50 "
   if ( (loopref_checked == "yes") or (looprem_checked == "yes") ):
      send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      #send_job_to_local_server(uniquedirname,rtcrcommand)
   else:
      #send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      send_job_to_local_server(uniquedirname,rtcrcommand)

   #find template for rendering in html
   #create tcr json file with template and sequence details
   a_segs = get_cdr_from_seq_by_aho_num(a_vd_seq, "A", True)
   aseq_fw = a_segs[5]
   aseq_cdr1 = a_segs[0]
   aseq_cdr2hv4 = a_segs[2]
   aseq_cdr3 = a_segs[4]
   aseq_cdr3_extnd = a_segs[6]
   tcrdata['aseq_vdomain'] = a_vd_seq
   tcrdata['aseq_fw'] = aseq_fw
   tcrdata['aseq_cdr1'] = aseq_cdr1
   tcrdata['aseq_cdr2hv4'] = aseq_cdr2hv4 
   tcrdata['aseq_cdr3'] = aseq_cdr3
   tcrdata['aseq_cdr3_extnd'] = aseq_cdr3_extnd
   tcrdata['a_fw1'] = a_segs[7]
   tcrdata['a_fw2'] = a_segs[8]
   tcrdata['a_fw3'] = a_segs[9]
   tcrdata['a_fw4'] = a_segs[10]

   b_segs = get_cdr_from_seq_by_aho_num(b_vd_seq, "B", True)
   bseq_fw = b_segs[5]
   bseq_cdr1 = b_segs[0]
   bseq_cdr2hv4 = b_segs[2]
   bseq_cdr3 = b_segs[4]
   bseq_cdr3_extnd = b_segs[6]
   tcrdata['bseq_vdomain'] = b_vd_seq
   tcrdata['bseq_fw'] = bseq_fw
   tcrdata['bseq_cdr1'] = bseq_cdr1
   tcrdata['bseq_cdr2hv4'] = bseq_cdr2hv4
   tcrdata['bseq_cdr3'] = bseq_cdr3
   tcrdata['bseq_cdr3_extnd'] = bseq_cdr3_extnd
   tcrdata['b_fw1'] = b_segs[7]
   tcrdata['b_fw2'] = b_segs[8]
   tcrdata['b_fw3'] = b_segs[9]
   tcrdata['b_fw4'] = b_segs[10]

   orientation_template_file = tmplt_db_seq+"TCR_FW_ORIENTATION.seq"
   ori_tmplt = find_orientation_template(aseq_fw.strip(), bseq_fw.strip(), orientation_template_file, None, pdb_blacklist)
   if all(ori_tmplt):
      tcrdata['ori_tmplt'] = str(ori_tmplt[0])[:4].upper()
      tcrdata['ori_tmplt_Apdb'] = str(ori_tmplt[0])[:6]
      tcrdata['ori_tmplt_Bpdb'] = str(ori_tmplt[1])[:6]
   else:
      errormsg = "No Template found for orientaion of Alpha and Beta chain sequences"
      return render_template("error.html", errormsg=errormsg)

   gm_multidb = glob.glob(tmplt_db_seq+"*_TCR_GM.fasta")
   vd_multidb = glob.glob(tmplt_db_seq+"*_TCR_VD.fasta")
   cdr1_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
   cdr2hv4_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
   agm_tmplt = seq_match_from_fasta_files_list( a_gm_seq, gm_multidb, pdb_blacklist )
   if agm_tmplt is not None:
      tcrdata['acdr1_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr1_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr1_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      tcrdata['acdr2hv4_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr2hv4_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr2hv4_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      #agm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), vd_multidb)
      #agm_tmplt_segs = get_cdr_from_seq_by_aho_num(agm_tmplt_vd_seq, "A", True)
      #tcrdata['acdr1_tmplt_seq'] = str(agm_tmplt_segs[0])
      #tcrdata['acdr2hv4_tmplt_seq'] = str(agm_tmplt_segs[2])
      agm_acdr1_tmplt_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), cdr1_multidb)
      tcrdata['acdr1_tmplt_seq'] = str(agm_acdr1_tmplt_seq)
      agm_acdr2hv4_tmplt_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), cdr2hv4_multidb)
      tcrdata['acdr2hv4_tmplt_seq'] = str(agm_acdr2hv4_tmplt_seq)
   else:
      multidb = glob.glob(tmplt_db_seq+"*_TCR_FW.fasta")
      find_templates(aseq_fw, multidb, "afw", "Alpha Framework", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
      find_templates(aseq_cdr1, multidb, "acdr1", "Alpha CDR1", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
      find_templates(aseq_cdr2hv4, multidb, "acdr2hv4", "Alpha CDR2 & HV4", None, pdb_blacklist)
   bgm_tmplt = seq_match_from_fasta_files_list( b_gm_seq, gm_multidb, pdb_blacklist )
   if bgm_tmplt is not None:
      tcrdata['bcdr1_tmplt_pdb'] = str(bgm_tmplt.id)[:6].upper()
      tcrdata['bcdr1_tmplt_pdbid'] = str(bgm_tmplt.id)[:4].upper()
      tcrdata['bcdr1_tmplt_pdb_chain'] = str(bgm_tmplt.id)[5:6].upper()
      tcrdata['bcdr2hv4_tmplt_pdb'] = str(bgm_tmplt.id)[:6].upper()
      tcrdata['bcdr2hv4_tmplt_pdbid'] = str(bgm_tmplt.id)[:4].upper()
      tcrdata['bcdr2hv4_tmplt_pdb_chain'] = str(bgm_tmplt.id)[5:6].upper()
      #bgm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), vd_multidb)
      #bgm_tmplt_segs = get_cdr_from_seq_by_aho_num(bgm_tmplt_vd_seq, "B", True)
      #tcrdata['bcdr1_tmplt_seq'] = str(bgm_tmplt_segs[0])
      #tcrdata['bcdr2hv4_tmplt_seq'] = str(bgm_tmplt_segs[2])
      bgm_acdr1_tmplt_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), cdr1_multidb)
      tcrdata['bcdr1_tmplt_seq'] = str(bgm_acdr1_tmplt_seq)
      bgm_acdr2hv4_tmplt_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), cdr2hv4_multidb)
      tcrdata['bcdr2hv4_tmplt_seq'] = str(bgm_acdr2hv4_tmplt_seq)
   else:
      multidb = glob.glob(tmplt_db_seq+"*_TCR_FW.fasta")
      find_templates(bseq_fw, multidb, "bfw", "Beta Framework", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
      find_templates(bseq_cdr1, multidb, "bcdr1", "Beta CDR1", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
      find_templates(bseq_cdr2hv4, multidb, "bcdr2hv4", "Beta CDR2 & HV4", None, pdb_blacklist)

   #CDR3
   cdr3_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR3_*.fasta")
   find_templates(aseq_cdr3, cdr3_multidb, "acdr3", "Alpha CDR3", None, pdb_blacklist)
   find_templates(bseq_cdr3, cdr3_multidb, "bcdr3", "Beta CDR3", None, pdb_blacklist)
   #override cdr3 seq with cdr3_extnd seq
   #cdr3a_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['acdr3_tmplt_id'], vd_multidb)
   #print  "cdr3a_tmplt_vd_seq", cdr3a_tmplt_vd_seq
   #cdr3a_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3a_tmplt_vd_seq, None, True)
   #tcrdata['acdr3_tmplt_seq'] = str(cdr3a_tmplt_segs[6])
   #cdr3b_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['bcdr3_tmplt_id'], vd_multidb)
   #cdr3b_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3b_tmplt_vd_seq, None, True)
   #tcrdata['bcdr3_tmplt_seq'] = str(cdr3b_tmplt_segs[6])
   cdr3extnd_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR3extnd_*.fasta")
   acdr3extnd_tmplt_seq = seq_by_id_from_fasta_files_list(tcrdata['acdr3_tmplt_id'], cdr3extnd_multidb)
   tcrdata['acdr3_tmplt_seq'] = str(acdr3extnd_tmplt_seq)
   bcdr3extnd_tmplt_seq = seq_by_id_from_fasta_files_list(tcrdata['bcdr3_tmplt_id'], cdr3extnd_multidb)
   tcrdata['bcdr3_tmplt_seq'] = str(bcdr3extnd_tmplt_seq)

   with open(uniquedirname+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")

   return redirect(url_for('rtcr', jobid=uniquedirname))


@app.route('/process_tcrpmhc_job/<aseq>/<bseq>/<pseq>/<mhc1aseq>/<mhc2aseq>/<mhc2bseq>/<loopref_checked>/<pdb_blacklist>')
def process_tcrpmhc_job(aseq,bseq,pseq,mhc1aseq,mhc2aseq,mhc2bseq,loopref_checked,pdb_blacklist):
   aseq=''.join(i for i in aseq if i.isalpha())
   bseq=''.join(i for i in bseq if i.isalpha())
   tcrdata['aseq_user']=aseq
   tcrdata['bseq_user']=bseq

   a_regexp_res = get_vdomain_using_regex(aseq,"A") 
   if (not a_regexp_res):
      errormsg = "No TCR Variable domain sequence identified from the Alpha chain input! If you used 'Select from gene', the CDR3 sequence was invalid."
      return render_template("error.html", errormsg=errormsg)
   a_vd_seq = a_regexp_res.groups()[0]
   a_gm_seq = a_regexp_res.groups()[1]+a_regexp_res.groups()[2]

   b_regexp_res = get_vdomain_using_regex(bseq,"B") 
   if (not b_regexp_res):
      errormsg = "No TCR variable domain sequence identified from the Beta chain input! If you used 'Select from gene', the CDR3 sequence was invalid."
      return render_template("error.html", errormsg=errormsg)
   b_vd_seq = b_regexp_res.groups()[0]
   b_gm_seq = b_regexp_res.groups()[1]+b_regexp_res.groups()[2]

   #we use uniquedirname to differentiate mhc1 and mhc2 complexes in res_tcrpmhc function
   if (not mhc1aseq == 'NA'):
      uniquedirname = create_and_cd_to_unique_dir("MHC1")#TCRC for tcr-p-mhc complex
      pseq=''.join(i for i in pseq if i.isalpha())
      mhc1aseq=''.join(i for i in mhc1aseq if i.isalpha())
      tcrdata['pep'] = pseq
      tcrdata['mhc1a'] = mhc1aseq
      rtcrcommand = "-blastp_path /www/cgi-bin/ncbi-blast-2.12.0+/bin/blastp -mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s -peptide %s -mhc1_alpha %s "  % (aseq, bseq, pseq, mhc1aseq) 
   elif ( (not mhc2aseq == 'NA') and (not mhc2bseq == 'NA') ):
      uniquedirname = create_and_cd_to_unique_dir("MHC2")#TCRC for tcr-p-mhc complex
      core_pseq = get_core_binding_pepseq(pseq)
      mhc2aseq=''.join(i for i in mhc2aseq if i.isalpha())
      mhc2bseq=''.join(i for i in mhc2bseq if i.isalpha())
      tcrdata['trunc_pep'] = core_pseq
      tcrdata['mhc2a'] = mhc2aseq
      tcrdata['mhc2b'] = mhc2bseq
      rtcrcommand = "-blastp_path /www/cgi-bin/ncbi-blast-2.12.0+/bin/blastp -mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s -peptide %s -mhc2_alpha %s -mhc2_beta %s "  % (aseq, bseq, core_pseq, mhc2aseq, mhc2bseq)  

   tcrdata['jobid'] = uniquedirname
   if pdb_blacklist:
      pdb_blacklist = pdb_blacklist.split(',')
      with open("ignore_list.txt", 'w') as ignore_list_file:
         for ignorepdb in pdb_blacklist:
            ignore_list_file.write(ignorepdb+'\n')
            rtcrcommand += "-ignore_list ignore_list.txt "
   #use -include_list include_list.txt flag to model on specific template
   #create a file with specific template name
   #rtcrcommand += "-include_list include_list.txt "

   looprem_checked = "no";
   if (looprem_checked == "yes"):
      rtcrcommand += "-remodel_tcr_cdr3_loops "
   if (loopref_checked == "yes"):
      rtcrcommand += "-refine_tcr_cdr3_loops -loops::max_inner_cycles 50 "
   if ( (loopref_checked == "yes") or (looprem_checked == "yes") ):
      #send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      send_tcrpmhc_job_to_local_server(uniquedirname,rtcrcommand)
   else:
      #send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      send_tcrpmhc_job_to_local_server(uniquedirname,rtcrcommand)

   #find template for rendering in html
   #create tcr json file with template and sequence details
   a_segs = get_cdr_from_seq_by_aho_num(a_vd_seq, "A", True)
   aseq_fw = a_segs[5]
   aseq_cdr1 = a_segs[0]
   aseq_cdr2hv4 = a_segs[2]
   aseq_cdr3 = a_segs[4]
   aseq_cdr3_extnd = a_segs[6]
   tcrdata['aseq_vdomain'] = a_vd_seq
   tcrdata['aseq_fw'] = aseq_fw
   tcrdata['aseq_cdr1'] = aseq_cdr1
   tcrdata['aseq_cdr2hv4'] = aseq_cdr2hv4 
   tcrdata['aseq_cdr3'] = aseq_cdr3
   tcrdata['aseq_cdr3_extnd'] = aseq_cdr3_extnd
   tcrdata['a_fw1'] = a_segs[7]
   tcrdata['a_fw2'] = a_segs[8]
   tcrdata['a_fw3'] = a_segs[9]
   tcrdata['a_fw4'] = a_segs[10]

   b_segs = get_cdr_from_seq_by_aho_num(b_vd_seq, "B", True)
   bseq_fw = b_segs[5]
   bseq_cdr1 = b_segs[0]
   bseq_cdr2hv4 = b_segs[2]
   bseq_cdr3 = b_segs[4]
   bseq_cdr3_extnd = b_segs[6]
   tcrdata['bseq_vdomain'] = b_vd_seq
   tcrdata['bseq_fw'] = bseq_fw
   tcrdata['bseq_cdr1'] = bseq_cdr1
   tcrdata['bseq_cdr2hv4'] = bseq_cdr2hv4
   tcrdata['bseq_cdr3'] = bseq_cdr3
   tcrdata['bseq_cdr3_extnd'] = bseq_cdr3_extnd
   tcrdata['b_fw1'] = b_segs[7]
   tcrdata['b_fw2'] = b_segs[8]
   tcrdata['b_fw3'] = b_segs[9]
   tcrdata['b_fw4'] = b_segs[10]

   orientation_template_file = tmplt_db_seq+"TCR_FW_ORIENTATION.seq"
   ori_tmplt = find_orientation_template(aseq_fw.strip(), bseq_fw.strip(), orientation_template_file, None, pdb_blacklist)
   if all(ori_tmplt):
      tcrdata['ori_tmplt'] = str(ori_tmplt[0])[:4].upper()
      tcrdata['ori_tmplt_Apdb'] = str(ori_tmplt[0])[:6]
      tcrdata['ori_tmplt_Bpdb'] = str(ori_tmplt[1])[:6]
   else:
      errormsg = "No Template found for orientaion of Alpha and Beta chain sequences"
      return render_template("error.html", errormsg=errormsg)
   gm_multidb = glob.glob(tmplt_db_seq+"*_TCR_GM.fasta")
   vd_multidb = glob.glob(tmplt_db_seq+"*_TCR_VD.fasta")
   cdr1_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
   cdr2hv4_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
   agm_tmplt = seq_match_from_fasta_files_list( a_gm_seq, gm_multidb, pdb_blacklist )
   if agm_tmplt is not None:
      tcrdata['acdr1_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr1_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr1_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      tcrdata['acdr2hv4_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr2hv4_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr2hv4_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      #agm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), vd_multidb)
      #agm_tmplt_segs = get_cdr_from_seq_by_aho_num(agm_tmplt_vd_seq, "A", True)
      #tcrdata['acdr1_tmplt_seq'] = str(agm_tmplt_segs[0])
      #tcrdata['acdr2hv4_tmplt_seq'] = str(agm_tmplt_segs[2])
      agm_acdr1_tmplt_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), cdr1_multidb)
      tcrdata['acdr1_tmplt_seq'] = str(agm_acdr1_tmplt_seq)
      agm_acdr2hv4_tmplt_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), cdr2hv4_multidb)
      tcrdata['acdr2hv4_tmplt_seq'] = str(agm_acdr2hv4_tmplt_seq)
   else:
      multidb = glob.glob(tmplt_db_seq+"*_TCR_FW.fasta")
      find_templates(aseq_fw, multidb, "afw", "Alpha Framework", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
      find_templates(aseq_cdr1, multidb, "acdr1", "Alpha CDR1", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
      find_templates(aseq_cdr2hv4, multidb, "acdr2hv4", "Alpha CDR2 & HV4", None, pdb_blacklist)
   bgm_tmplt = seq_match_from_fasta_files_list( b_gm_seq, gm_multidb, pdb_blacklist )
   if bgm_tmplt is not None:
      tcrdata['bcdr1_tmplt_pdb'] = str(bgm_tmplt.id)[:6].upper()
      tcrdata['bcdr1_tmplt_pdbid'] = str(bgm_tmplt.id)[:4].upper()
      tcrdata['bcdr1_tmplt_pdb_chain'] = str(bgm_tmplt.id)[5:6].upper()
      tcrdata['bcdr2hv4_tmplt_pdb'] = str(bgm_tmplt.id)[:6].upper()
      tcrdata['bcdr2hv4_tmplt_pdbid'] = str(bgm_tmplt.id)[:4].upper()
      tcrdata['bcdr2hv4_tmplt_pdb_chain'] = str(bgm_tmplt.id)[5:6].upper()
      #bgm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), vd_multidb)
      #bgm_tmplt_segs = get_cdr_from_seq_by_aho_num(bgm_tmplt_vd_seq, "B", True)
      #tcrdata['bcdr1_tmplt_seq'] = str(bgm_tmplt_segs[0])
      #tcrdata['bcdr2hv4_tmplt_seq'] = str(bgm_tmplt_segs[2])
      bgm_acdr1_tmplt_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), cdr1_multidb)
      tcrdata['bcdr1_tmplt_seq'] = str(bgm_acdr1_tmplt_seq)
      bgm_acdr2hv4_tmplt_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), cdr2hv4_multidb)
      tcrdata['bcdr2hv4_tmplt_seq'] = str(bgm_acdr2hv4_tmplt_seq)
   else:
      multidb = glob.glob(tmplt_db_seq+"*_TCR_FW.fasta")
      find_templates(bseq_fw, multidb, "bfw", "Beta Framework", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
      find_templates(bseq_cdr1, multidb, "bcdr1", "Beta CDR1", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
      find_templates(bseq_cdr2hv4, multidb, "bcdr2hv4", "Beta CDR2 & HV4", None, pdb_blacklist)

   #CDR3
   cdr3_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR3_*.fasta")
   find_templates(aseq_cdr3, cdr3_multidb, "acdr3", "Alpha CDR3", None, pdb_blacklist)
   find_templates(bseq_cdr3, cdr3_multidb, "bcdr3", "Beta CDR3", None, pdb_blacklist)
   #override cdr3 seq with cdr3_extnd seq

   #cdr3a_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['acdr3_tmplt_id'], vd_multidb)
   #print  "cdr3a_tmplt_vd_seq", cdr3a_tmplt_vd_seq
   #cdr3a_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3a_tmplt_vd_seq, None, True)
   #tcrdata['acdr3_tmplt_seq'] = str(cdr3a_tmplt_segs[6])
   #cdr3b_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['bcdr3_tmplt_id'], vd_multidb)
   #cdr3b_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3b_tmplt_vd_seq, None, True)
   #tcrdata['bcdr3_tmplt_seq'] = str(cdr3b_tmplt_segs[6])
   cdr3extnd_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR3extnd_*.fasta")
   acdr3extnd_tmplt_seq = seq_by_id_from_fasta_files_list(tcrdata['acdr3_tmplt_id'], cdr3extnd_multidb)
   tcrdata['acdr3_tmplt_seq'] = str(acdr3extnd_tmplt_seq)
   bcdr3extnd_tmplt_seq = seq_by_id_from_fasta_files_list(tcrdata['bcdr3_tmplt_id'], cdr3extnd_multidb)
   tcrdata['bcdr3_tmplt_seq'] = str(bcdr3extnd_tmplt_seq)

   with open(uniquedirname+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")
   return redirect(url_for('res_tcrpmhc', jobid=uniquedirname))

def send_job_to_local_server(uniquedirname,rtcrcommand):
   rtcrcommand +=   " -tcr_template_db_path /www/tcrmodel/static/downloads/additional_protocol_data/tcr/"
   commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/tcrmodel.static.linuxgccrelease -database /www/cgi-bin/rosetta/Rosetta/main/database " + rtcrcommand +  " > res.out "
   p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   #out,err = p.communicate()
   #Do not call communicate if you do not want to wait for the command to complete
   return 

def send_tcrpmhc_job_to_local_server(uniquedirname,rtcrcommand):
   rtcrcommand +=   " -tcr_template_db_path /www/tcrmodel/static/downloads/additional_protocol_data/tcr/"
   commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/tcr_complex_model.static.linuxgccrelease -database /www/cgi-bin/rosetta/Rosetta/main/database " + rtcrcommand +  " > res.out "
   commandline2 = "echo " + commandline + " > ros_cmd.txt"
   p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   p2 = subprocess.Popen(commandline2, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   #out,err = p.communicate()
   #Do not call communicate if you do not want to wait for the command to complete
   #print stdout
   #print stderr
   return 

def send_job_to_ibbr_cluster(uniquejobid,rtcrcommand):
   jobscriptfile = "RosettaTCR.slurm"
   pbs_file_path = os.path.join(rundir_path,uniquejobid,jobscriptfile)
   ig_fn = "ignore_list.txt"
   with open(pbs_file_path, 'a') as f:
      f.write('#!/bin/sh\n')
      f.write('#SBATCH -n 16\n')
      f.write('#SBATCH --time 72:00:00\n')
      f.write('#SBATCH --mem=512M\n')
      f.write('cd $SLURM_SUBMIT_DIR\n')
      f.write('rundir=${PWD##*/}\n')
      commandline = "/home/gowthamanr/Rosetta/main/source/bin/tcrmodel.static.linuxgccrelease " + rtcrcommand +  " > res.out "
      f.write(commandline+'\n')
      f.write("if test -f res.out ; then scp res.out gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/res.out ;fi\n")
      f.write("chmod 777 $SLURM_SUBMIT_DIR/tcrmodel.pdb\n")
      f.write("if test -f tcrmodel.pdb ; then scp tcrmodel.pdb gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcrmodel.pdb ;fi\n")
      f.write("if ! test -f tcrmodel.pdb ; then echo ERROR > tcr.fail; scp tcr.fail gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcr.fail ;fi\n")
   submit_job_on_cluster(uniquejobid,jobscriptfile,ig_fn)
   return 

def renumber_tcrpdb_by_aho_number(outpdb):
   outpdb_orig = 'tcrmodel.pdb.orig'
   tcr_achain = 'tcr_achain.pdb'
   tcr_bchain = 'tcr_bchain.pdb'
   tcr_achain_aho = 'tcr_achain_aho.pdb'
   tcr_bchain_aho = 'tcr_bchain_aho.pdb'
   make_pdb(outpdb, "A", tcr_achain)
   make_pdb(outpdb, "B", tcr_bchain)
   renumber_pdbfile_to_aho(tcr_achain, "A", tcr_achain_aho, True)
   renumber_pdbfile_to_aho(tcr_bchain, "B", tcr_bchain_aho, True)
   if ( (os.path.isfile(tcr_achain_aho) and os.path.getsize(tcr_achain_aho) > 0)
        and (os.path.isfile(tcr_bchain_aho) and os.path.getsize(tcr_bchain_aho) > 0) ):
      shutil.move(outpdb, outpdb_orig)
      with open(outpdb,'wb') as wfd:
         for f in [tcr_achain_aho, tcr_bchain_aho]:
            with open(f,'rb') as fd:
               shutil.copyfileobj(fd, wfd, 1024*1024*10)
   return

def add_ss_header_to_pdbfile(outpdb):
   dsspfile = "tcrmodel.dssp"
   ssheaderfile = "tcrmodel_ssheader.pdb"
   tcr_ss_tmp_pdb = "tcrmodel.tmp.pdb"
   dssp2pdb_script = os.path.join(app.root_path, "dssp2pdb.pl")
   subprocess.call(["/www/cgi-bin/dssp-3.0.0/mkdssp", outpdb, dsspfile])
   if ( os.path.isfile(dsspfile) and (os.path.getsize(dsspfile)>0) ):
      with open(ssheaderfile, 'w') as outfile:
         subprocess.call(["perl", dssp2pdb_script, "-35", dsspfile], stdout=outfile)
   if ( os.path.isfile(ssheaderfile) and (os.path.getsize(ssheaderfile)>0) ):
      with open(tcr_ss_tmp_pdb,'w') as wfd:
         for f in [ssheaderfile, outpdb]:
            with open(f,'rb') as fd:
               shutil.copyfileobj(fd, wfd)
      shutil.move(tcr_ss_tmp_pdb, outpdb)         
   return

def create_profit_infile_cdr1(infile,refchain,mobchain,outpdb):
   f = open(infile,"w")
   f.write("ATOMS N,CA,C\n")
   f.write("ZONE CLEAR\n")
   f.write("ZONE "+refchain+"23"+"-"+refchain+"44"+":"+mobchain+"23"+"-"+mobchain+"44"+"\n")
   f.write("FIT\n")
   f.write("WRITE "+outpdb+"\n")
   f.close()

def create_profit_infile_cdr3(infile,refchain,mobchain,outpdb):
   f = open(infile,"w")
   f.write("ATOMS N,CA,C\n")
   f.write("ZONE CLEAR\n")
   f.write("ZONE "+refchain+"107"+"-"+refchain+"138"+":"+mobchain+"107"+"-"+mobchain+"138"+"\n")
   f.write("FIT\n")
   f.write("WRITE "+outpdb+"\n")
   f.close()

def create_profit_infile_cdr2hv4(infile,refchain,mobchain,outpdb):
   f = open(infile,"w")
   f.write("ATOMS N,CA,C\n")
   f.write("ZONE CLEAR\n")
   f.write("ZONE "+refchain+"56"+"-"+refchain+"70"+":"+mobchain+"56"+"-"+mobchain+"70"+"\n")
   f.write("FIT\n")
   f.write("WRITE "+outpdb+"\n")
   f.close()

def create_profit_infile_ori(infile,refchain,mobchain,outpdb):
   f = open(infile,"w")
   f.write("ATOMS N,CA,C\n")
   f.write("ZONE CLEAR\n")
   f.write("ZONE "+refchain+"4"+"-"+refchain+"23"+":"+mobchain+"4"+"-"+mobchain+"23"+"\n")
   f.write("ZONE "+refchain+"44"+"-"+refchain+"56"+":"+mobchain+"44"+"-"+mobchain+"56"+"\n")
   f.write("ZONE "+refchain+"70"+"-"+refchain+"107"+":"+mobchain+"70"+"-"+mobchain+"107"+"\n")
   f.write("ZONE "+refchain+"138"+"-"+refchain+"148"+":"+mobchain+"138"+"-"+mobchain+"148"+"\n")
   f.write("FIT\n")
   f.write("WRITE "+outpdb+"\n")
   f.close()

def align_tmplt_and_target(infile, reference, mobile):
   profit_program = "profit"
   cmd = subprocess.Popen([profit_program, '-f', infile, '-h', reference, mobile], stdout=subprocess.PIPE)
   return

def get_pdb_templates(tcrjsondata):
   ref = "tcrmodel.pdb"
   refchainA = "A"
   refchainB = "B"
   #cdr1a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr1_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr1_tmplt_pdb'][5:6]
   acdr1_outpdb = tcrjsondata['acdr1_tmplt_pdb'] + "_Acdr1_tmplt.pdb"
   infile = "Acdr1_profit.in"
   create_profit_infile_cdr1(infile,refchainA,mobchain,acdr1_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr2a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr2hv4_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr2hv4_tmplt_pdb'][5:6]
   acdr2hv4_outpdb = tcrjsondata['acdr2hv4_tmplt_pdb'] + "_Acdr2hv4_tmplt.pdb"
   infile = "Acdr2hv4_profit.in"
   create_profit_infile_cdr2hv4(infile,refchainA,mobchain,acdr2hv4_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr3a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr3_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr3_tmplt_pdb'][5:6]
   acdr3_outpdb = tcrjsondata['acdr3_tmplt_pdb'] + "_Acdr3_tmplt.pdb"
   infile = "Acdr3_profit.in"
   create_profit_infile_cdr3(infile,refchainA,mobchain,acdr3_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr1b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr1_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr1_tmplt_pdb'][5:6]
   bcdr1_outpdb = tcrjsondata['bcdr1_tmplt_pdb'] + "_Bcdr1_tmplt.pdb"
   infile = "Bcdr1_profit.in"
   create_profit_infile_cdr1(infile,refchainB,mobchain,bcdr1_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr2b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr2hv4_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr2hv4_tmplt_pdb'][5:6]
   bcdr2hv4_outpdb = tcrjsondata['bcdr2hv4_tmplt_pdb'] + "_Bcdr2hv4_tmplt.pdb"
   infile = "Bcdr2hv4_profit.in"
   create_profit_infile_cdr2hv4(infile,refchainB,mobchain,bcdr2hv4_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr3b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr3_tmplt_pdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr3_tmplt_pdb'][5:6]
   bcdr3_outpdb = tcrjsondata['bcdr3_tmplt_pdb'] + "_Bcdr3_tmplt.pdb"
   infile = "Bcdr3_profit.in"
   create_profit_infile_cdr3(infile,refchainB,mobchain,bcdr3_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #orientation templates
   #orientation A
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['ori_tmplt_Apdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['ori_tmplt_Apdb'][5:6]
   oria_outpdb = tcrjsondata['ori_tmplt_Apdb'] + "_oriA_tmplt.pdb"
   infile = "oriA_profit.in"
   create_profit_infile_ori(infile,refchainA,mobchain,oria_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #orientation B
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['ori_tmplt_Bpdb'][:4].lower()+"_aho.pdb.gz")[0]
   mobchain = tcrjsondata['ori_tmplt_Bpdb'][5:6]
   orib_outpdb = tcrjsondata['ori_tmplt_Bpdb'] + "_oriB_tmplt.pdb"
   infile = "oriB_profit.in"
   create_profit_infile_ori(infile,refchainB,mobchain,orib_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   '''   
   tar_filename = tcrjsondata['jobid'] + ".tar.gz"
   from contextlib import closing
   with closing(tarfile.open(tar_filename, "w:gz")) as tar:
   for name in [ "tcrmodel.pdb", acdr1_outpdb, acdr2hv4_outpdb, acdr3_outpdb, bcdr1_outpdb, bcdr2hv4_outpdb, bcdr3_outpdb, oria_outpdb,orib_outpdb]:
   tar.add(name)
   '''      
   return

@app.route('/rtcr/<jobid>')
def rtcr(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   os.chdir(outdir)
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   tcrjsonfile_spath = os.path.join(rundir_spath,str(jobid),str(jobid)+'.json')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   outpdb = os.path.join(rundir_path,str(jobid),'tcrmodel.pdb')
   errorfile = os.path.join(rundir_path,str(jobid),'tcr.fail')
   errorfile2 = os.path.join(rundir_path,str(jobid),'ROSETTA_CRASH.log')
   if os.path.isfile(outpdb):
      renumber_tcrpdb_by_aho_number(outpdb)
      add_ss_header_to_pdbfile(outpdb)
      outpdb_spath = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
      stdout = get_pdb_templates(tcrjsondata)
      return render_template("tcrmodel_completed.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   
   elif os.path.isfile(errorfile):
      return render_template("rtcrfailed.html", rtcrjobid = jobid)   
   elif os.path.isfile(errorfile2):
      return render_template("rtcrfailed.html", rtcrjobid = jobid)   
   else:
      return render_template("rtcrqueued.html", jobid=jobid, tj=tcrjsondata)   

@app.route('/res_tcrpmhc/<jobid>')
def res_tcrpmhc(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   os.chdir(outdir)
   modelpdb = os.path.join(rundir_path,str(jobid),'tcrpmhc_model.pdb')
   errorfile = os.path.join(rundir_path,str(jobid),'tcr.fail')
   errorfile2 = os.path.join(rundir_path,str(jobid),'ROSETTA_CRASH.log')
   if os.path.isfile(modelpdb):
      #run model delsection here?
      renamed_chain_file = "renamed_chain.pdb"
      renum_aho_file = "renum_aho.pdb"
      ss_header_file = "tcrmodel_ssheader.pdb"
      if not os.path.isfile(renamed_chain_file): 
         if (jobid[:4] == "MHC1"):
            rename_pdb_chains("ABCD", "ACDE", modelpdb, renamed_chain_file )
         if (jobid[:4] == "MHC2"):
            rename_pdb_chains("ABCDE", "ABCDE", modelpdb, renamed_chain_file )
      if not os.path.isfile(renum_aho_file): 
         renumber_pdbchain_to_aho(renamed_chain_file, "D", None, "ahoD.pdb", True)
         renumber_pdbchain_to_aho("ahoD.pdb", "E", None, renum_aho_file, True)
         if (os.path.isfile(renum_aho_file) and os.path.getsize(renum_aho_file) > 0):
            shutil.copyfile(renum_aho_file, modelpdb)
      if not os.path.isfile(ss_header_file):
         add_ss_header_to_pdbfile(modelpdb)
      outpdb_spath = os.path.join(rundir_spath,str(jobid),'tcrpmhc_model.pdb')
      #stdout = get_pdb_templates(tcrjsondata)
      #Read Json file and change empty '' values to 'NA' values
      tcrjsonfile = os.path.join(rundir_path,str(jobid),"tcr_template_info.json")
      with open(tcrjsonfile) as json_data:
         tcrjsondata = json.load(json_data)
         for key in tcrjsondata:
            if not tcrjsondata[key]:
               tcrjsondata[key] = 'NA';

      from datetime import datetime
      now = datetime.now()
      current_time = now.strftime("%H:%M:%S")
      print("Current Time =", current_time)
      return render_template("tcrpmhc_completed.html", jobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   
   elif os.path.isfile(errorfile):
      return render_template("error.html", jobid=jobid, errormsg="Modeling could not complete!")   
   elif os.path.isfile(errorfile2):
      errormsg = "Modeling could not complete!"
      with open(errorfile2) as fh:
         for line in fh:
            if line.startswith("ERROR:"):
               errormsg = line[7:] 
               break 
      return render_template("error.html", jobid=jobid, errormsg=errormsg)   
   else:
      return render_template("tcrpmhc_queued.html", jobid=jobid)

@app.route('/mtcr/<jobid>')
def mtcr(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   os.chdir(outdir)
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   subdirs = next(os.walk('.'))[1]

   alljobsdone = True
   i = 0
   for subdir in subdirs:
      outpdb = os.path.join(rundir_path,str(jobid),subdir,str(subdir)+"_tcrmodel.pdb")
      outfail = os.path.join(rundir_path,str(jobid),subdir,"tcr.fail")
      if os.path.isfile(outpdb):
         currstatus = "Completed"
         clean_outpdb = os.path.join(rundir_path,str(jobid),subdir,str(subdir)+"_tcrmodel.clean.pdb")
         if not os.path.isfile(clean_outpdb):
            #change to subdir and add aho num and secondary structure info
            os.chdir(os.path.join(rundir_path,str(jobid),subdir))
            shutil.copyfile(outpdb, clean_outpdb)
            print "aho", outpdb
            renumber_tcrpdb_by_aho_number(clean_outpdb)
            print "ss", outpdb
            add_ss_header_to_pdbfile(clean_outpdb)
            os.remove(outpdb)         
            shutil.copyfile(clean_outpdb, outpdb)
            os.chdir(outdir)
      elif os.path.isfile(outfail):
         currstatus = "Failed"
      else:
         currstatus = "Running"
         alljobsdone = False
      tcrjsondata['tcrmodel'][i]['status'] = currstatus
      i += 1

   if alljobsdone:
      zipf = zipfile.ZipFile(str(jobid)+'.zip', 'w', zipfile.ZIP_DEFLATED)#zip files
      pdblist = glob.glob("*/*_tcrmodel.pdb")
      for pdbfile in pdblist:
         outpdb = os.path.join(rundir_path,str(jobid),pdbfile)
         #renumber_tcrpdb_by_aho_number(outpdb)
         #add_ss_header_to_pdbfile(outpdb)
         zipf.write(pdbfile)
      zipf.close()

      return render_template("mtcrcompleted.html", mtcrjobid = jobid, tj=tcrjsondata,rundir_spath=rundir_spath)
   else:
      return render_template("mtcrrunning.html", mtcrjobid = jobid, tj=tcrjsondata,rundir_spath=rundir_spath)

@app.route('/mtcrex/<jobid>')
def mtcrex(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   os.chdir(outdir)
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   return render_template("mtcrcompleted.html", mtcrjobid = jobid, tj=tcrjsondata,rundir_spath=rundir_spath)

@app.route('/rtcrex/<jobid>')
def rtcrex(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   outpdb_spath = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   return render_template("tcrmodel_completed.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   

@app.route('/viewmodel/<jobid>')
def viewmodel(jobid):
   modelfname = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   return render_template("view.html", jobid=jobid, modelfname=modelfname, tj=tcrjsondata)

@app.route('/mviewmodel/<jobid>/<prefixtag>')
def mviewmodel(jobid,prefixtag):
   modelfname = os.path.join(rundir_spath,str(jobid),str(prefixtag),str(prefixtag)+'_tcrmodel.pdb')
   return render_template("mview.html", jobid=jobid, modelfname=modelfname)
   #return render_template("error.html", errormsg=modelfname)

@app.route('/searchid', methods=['POST', 'GET'])
def searchid():
   if not request.form.get("jobidinput"):
      return render_template("error.html", errormsg="Job ID not entered")
   jobid = request.form['jobidinput'].strip()
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   if (jobid[:4] == "TCRM"):
      os.chdir(outdir)
      pdblist = glob.glob("*/*_tcrmodel.pdb")
      if pdblist:
         return redirect(url_for('mtcrex', jobid=jobid))
      else:
         return redirect(url_for('mtcr', jobid=jobid))
   elif (jobid[:4] == "TCRS"):
      outpdb = os.path.join(rundir_path,str(jobid),'tcrmodel.pdb')
      if os.path.isfile(outpdb):
         return redirect(url_for('rtcrex', jobid=jobid))
      else:
         return redirect(url_for('rtcr', jobid=jobid))
   elif (jobid[:3] == "MHC"):
      return redirect(url_for('res_tcrpmhc', jobid=jobid))
   else:
      return render_template("error.html", errormsg="Job ID not found")

@app.route('/tcrmodel_submitjob1', methods = ['POST', 'GET'])
def tcrmodel_submitjob1():
   aseq = request.form['alphachain']
   bseq = request.form['betachain']
   if not aseq:
      return render_template("error.html", errormsg="TCR alpha sequence not entered!")
   if not bseq:
      return render_template("error.html", errormsg="TCR beta sequence not entered!") 
 
   #simil_cutoff = request.form.get('simcutoff')
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   #looprem_checked = "yes" if (request.form.get("lr"))else "no"
   #uniquedirname = processjob(aseq,bseq,simil_cutoff,loopref_checked,looprem_checked,pdb_blacklist)
   return redirect(url_for('processjob', aseq=aseq,bseq=bseq,loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

@app.route('/tcrmodel_submitjob2', methods=['POST', 'GET'])
def tcrmodel_submitjob2():
   aspecies = request.form.get('sele_aspecies')
   trav = request.form.get('sele_trav')
   traj = request.form.get('sele_traj')
   acdr = request.form['sele_acdr']
   bspecies = request.form.get('sele_bspecies')   
   trbv = request.form.get('sele_trbv')
   trbj = request.form.get('sele_trbj')
   bcdr = request.form['sele_bcdr']
   if not aspecies:
      return render_template("error.html", errormsg="Species for alpha chain not selected!")
   if not trav:
      return render_template("error.html", errormsg="TRAV gene not selected!")
   if not traj:
      return render_template("error.html", errormsg="TRAJ gene not selected!")
   if not acdr:
      return render_template("error.html", errormsg="CDR3 sequence for alpha chain not entered!")
   if not bspecies:
      return render_template("error.html", errormsg="Species for beta chain not selected!")
   if not trbv:
      return render_template("error.html", errormsg="TRBV gene not selected!")
   if not trbj:
      return render_template("error.html", errormsg="TRBJ gene not selected!")
   if not bcdr:
      return render_template("error.html", errormsg="CDR3 sequence for beta chain not entered!")
   
   aseq =  trav+acdr.upper()+traj
   bseq =  trbv+bcdr.upper()+trbj
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   return redirect(url_for('processjob', aseq=aseq,bseq=bseq,loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

@app.route('/tcrpmhc_submitjob1', methods=['POST', 'GET'])
def tcrpmhc_submitjob1():
   from datetime import datetime
   now = datetime.now()
   current_time = now.strftime("%H:%M:%S")
   print("Current Time =", current_time)
   aseq = request.form['alphachain']
   bseq = request.form['betachain']
   if not aseq:
      return render_template("error.html", errormsg="TCR alpha sequence not entered!")
   if not bseq:
      return render_template("error.html", errormsg="TCR beta sequence not entered!")
   
   pseq = request.form['pepchain']
   mhc1aseq = request.form['mhc1aseq']
   mhc2aseq = request.form['mhc2aseq']
   mhc2bseq = request.form['mhc2bseq']
   if not pseq:
      return render_template("error.html", errormsg="Peptide sequence not entered!")
   if not mhc1aseq and not mhc2aseq and not mhc2bseq:
      return render_template("error.html", errormsg="Enter either MHC I or MHC II sequence(s)!")
   if mhc1aseq and len(pseq) < 8:
      return render_template("error.html", errormsg="Currently MHC I complex modeling supports only 8-mer peptides!")
   if (mhc2aseq and not mhc2bseq) or (mhc2bseq and not mhc2aseq):   
      return render_template("error.html", errormsg="Enter both alpha and beta sequences for MHC II!")
   if mhc2aseq and mhc2bseq and len(pseq) < 9:
      return render_template("error.html", errormsg="Currently MHC II complex modeling supports only 9-mer peptides!")

   if not mhc1aseq:
      mhc1aseq = 'NA'
   if not mhc2aseq:
      mhc2aseq = 'NA'
   if not mhc2bseq:
      mhc2bseq = 'NA'

   fields = [pseq, mhc1aseq, mhc2aseq, mhc2bseq]
   seq_names = ["Peptide", "MHC I", "MHC II alpha", "MHC II beta"]
   for i in range(len(fields)):
      if not re.search("^[ARNDCEQGHILKMFPSTWYV]+$", fields[i]):
         return render_template("error.html", errormsg=seq_names[i]+" sequence can only contain standard amino acids!") 

   #simil_cutoff = request.form.get('simcutoff')
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   #looprem_checked = "yes" if (request.form.get("lr"))else "no"
   #uniquedirname = processjob(aseq,bseq,simil_cutoff,loopref_checked,looprem_checked,pdb_blacklist)
   return redirect(url_for('process_tcrpmhc_job', aseq=aseq, bseq=bseq, pseq=pseq, mhc1aseq=mhc1aseq, mhc2aseq=mhc2aseq, mhc2bseq=mhc2bseq, loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

@app.route('/tcrpmhc_submitjob2', methods=['POST', 'GET'])
def tcrpmhc_submitjob2():
   aspecies = request.form.get('aspecies')
   trav = request.form.get('trav')
   traj = request.form.get('traj')
   acdr = request.form.get('acdr')
   bspecies = request.form.get('bspecies')
   trbv = request.form.get('trbv')
   trbj = request.form.get('trbj')
   bcdr = request.form.get('bcdr')
   if not aspecies:
      return render_template("error.html", errormsg="Species for alpha chain not selected!")
   if not trav:
      return render_template("error.html", errormsg="TRAV gene not selected!")
   if not traj:
      return render_template("error.html", errormsg="TRAJ gene not selected!")
   if not acdr:
      return render_template("error.html", errormsg="CDR3 sequence for alpha chain not entered!")
   if not bspecies:
      return render_template("error.html", errormsg="Species for beta chain not selected!")
   if not trbv:
      return render_template("error.html", errormsg="TRBV gene not selected!")
   if not trbj:
      return render_template("error.html", errormsg="TRBJ gene not selected!")
   if not bcdr:
      return render_template("error.html", errormsg="CDR3 sequence for beta chain not entered!")

   pseq = request.form['pepchain']
   mhc1aseq = request.form['mhc1aseq']
   mhc2aseq = request.form['mhc2aseq']
   mhc2bseq = request.form['mhc2bseq']
   if not pseq:
      return render_template("error.html", errormsg="Peptide sequence not entered!")
   if not mhc1aseq and not mhc2aseq and not mhc2bseq:
      return render_template("error.html", errormsg="Enter either MHC I or MHC II sequence(s)!")
   if mhc1aseq and len(pseq) < 8:
      return render_template("error.html", errormsg="Currently MHC I complex modeling supports only 8-mer peptides!")
   if (mhc2aseq and not mhc2bseq) or (mhc2bseq and not mhc2aseq):   
      return render_template("error.html", errormsg="Enter both alpha and beta sequences for MHC II!")
   if mhc2aseq and mhc2bseq and len(pseq) < 9:
      return render_template("error.html", errormsg="Currently MHC II complex modeling supports only 9-mer peptides!")
   
   if not mhc1aseq:
      mhc1aseq = 'NA'
   if not mhc2aseq:
      mhc2aseq = 'NA'
   if not mhc2bseq:
      mhc2bseq = 'NA'

   fields = [pseq, mhc1aseq, mhc2aseq, mhc2bseq]
   seq_names = ["Peptide", "MHC I", "MHC II alpha", "MHC II beta"]
   for i in range(len(fields)):
      if not re.search("^[ARNDCEQGHILKMFPSTWYV]+$", fields[i]):
         return render_template("error.html", errormsg=seq_names[i]+" sequence can only contain standard amino acids!") 
   
   aseq =  trav+acdr.upper()+traj
   bseq =  trbv+bcdr.upper()+trbj
   #simil_cutoff = request.form.get('simcutoff')
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   #looprem_checked = "yes" if (request.form.get("lr"))else "no"
   #uniquedirname = processjob(aseq,bseq,simil_cutoff,loopref_checked,looprem_checked,pdb_blacklist)
   return redirect(url_for('process_tcrpmhc_job', aseq=aseq, bseq=bseq, pseq=pseq, mhc1aseq=mhc1aseq, mhc2aseq=mhc2aseq, mhc2bseq=mhc2bseq, loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

@app.route("/DownloadFile/<path:filepath>")
def DownloadFile (filepath = None):
   if filepath is None:
      self.Error(400)
   try:
      return send_file(filepath, as_attachment=True)
   except Exception as e:
      self.log.exception(e)
      self.Error(400)

@app.route('/GetFilesFromFolder/<jobid>/<filename>')
def GetFilesFromFolder(jobid,filename):
   filepath = os.path.join(rundir_path,jobid)
   return send_from_directory(filepath, filename, as_attachment=True)

def allowed_file(filename):
   ALLOWED_EXTENSIONS = set(['txt', 'fasta', 'fa',])
   return '.' in filename and \
       filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def send_batchjob_to_ibbr_cluster(afile,bfile,loopref_checked,pdb_blacklist,uniquejobid):
   looprem_checked = "no";
   masterscript = "RosettaTCR.sh"
   ms_file_path = os.path.join(rundir_path,uniquejobid,masterscript)
   ig_fn = "ignore_list.txt"
   ig_file_path = os.path.join(rundir_path,uniquejobid,ig_fn)

   if pdb_blacklist:
      pdb_blacklist = pdb_blacklist.split(',')
      with open("ignore_list.txt", 'w') as ignore_list_file:
         for ignorepdb in pdb_blacklist:
            ignore_list_file.write(ignorepdb+'\n')
            
   arecords = list(SeqIO.parse(afile, "fasta"))
   brecords = list(SeqIO.parse(bfile, "fasta"))
   tcrdata['tcrmodel'] = []  
         
   with open(ms_file_path, 'a') as ms_file:
      ms_file.write('#!/bin/sh\n')
      pbs_err_file = uniquejobid + ".err"
      pbs_out_file = uniquejobid + ".out"
      #pbs_file.write('cd $PBS_O_WORKDIR\n')
      #pbs_file.write('rundir=${PWD##*/}\n')
      prefixnum = 0
      tcrdata['jobid'] = uniquejobid
      for (arecord,brecord) in zip(arecords, brecords):
         #check if the input tcr sequence can be parsed, otherwise skip the sequence
         aregexres = get_vdomain_using_regex(arecord.seq, "A")
         bregexres = get_vdomain_using_regex(brecord.seq, "B")
         if not aregexres: continue
         if not bregexres: continue
         prefixnum += 1
         prefixtag = str(arecord.name[:4].strip()) + "_" + str(prefixnum)
         jobdir = os.path.join(rundir_path,uniquejobid,prefixtag)
         if not os.path.exists(jobdir):
            os.makedirs(jobdir)
            os.chmod(jobdir, 0777)
         jobfile = prefixtag+".sh"
         jobfilepath = os.path.join(rundir_path,uniquejobid,prefixtag,jobfile)
         f = open(jobfilepath, "w")
         ms_file.write('cd '+prefixtag+'\n') 
         ms_file.write('sbatch '+jobfile+'\n') 
         ms_file.write('cd ..\n') 
         aregexres = get_vdomain_using_regex(arecord.seq, "A")
         bregexres = get_vdomain_using_regex(brecord.seq, "B")
         tcrdata['tcrmodel'].append({  
               'prefixnum': prefixnum,
               'prefixtag': prefixtag,
               'aseq_user': str(arecord.seq.strip()),
               'bseq_user': str(brecord.seq.strip()),
               'aseq_vdomain': aregexres.groups()[0],
               'bseq_vdomain': bregexres.groups()[0],
               'status': ''
         })
         rtcrcommand = "-mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s "  % (arecord.seq,brecord.seq) 
         #rtcrcommand = "-blastp_identity_cutoff 90 -mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s "  % (arecord.seq,brecord.seq) 
         rtcrcommand += " -out:prefix " + prefixtag + "_ "
         if (loopref_checked == "yes"):
            rtcrcommand += "-refine_tcr_cdr3_loops "
         if pdb_blacklist:
            rtcrcommand += "-ignore_list ignore_list.txt "
         commandline = "/home/gowthamanr/Rosetta/main/source/bin/tcrmodel.static.linuxgccrelease " + rtcrcommand +  " > res.out "
         f.write("#!/bin/sh\n")
         f.write('#SBATCH -n 16\n')
         f.write('#SBATCH --time 72:00:00\n')
         f.write('#SBATCH --mem=512M\n')
         f.write('cd $SLURM_SUBMIT_DIR\n')
         f.write('rundir=${PWD##*/}\n')
         f.write(commandline+'\n')
         res_file = "res.out"
         out_pdb = prefixtag + "_tcrmodel.pdb"
         out_fail = "tcr.fail"
         out_res_file = os.path.join(rundir_path,uniquejobid,prefixtag,res_file)
         out_pdb_file = os.path.join(rundir_path,uniquejobid,prefixtag,out_pdb)
         out_fail_file = os.path.join(rundir_path,uniquejobid,prefixtag,out_fail)
         f.write("chmod 777 "+res_file+'\n')
         f.write("if test -f " + res_file + " ; then scp " + res_file + " tcrmodel@piercelab.ibbr.umd.edu:"+out_res_file+" ;fi\n")
         f.write("if test -f " + out_pdb + " ; then chmod 755 " + out_pdb + "; scp " + out_pdb + " tcrmodel@piercelab.ibbr.umd.edu:"+out_pdb_file+" ;fi\n")
         f.write("if ! test -f " + out_pdb + " ; then echo ERROR > tcr.fail; chmod 755 tcr.fail; scp tcr.fail tcrmodel@piercelab.ibbr.umd.edu:"+out_fail_file+" ;fi\n")
         f.close()
   submit_batchjob_on_cluster(uniquejobid,masterscript,ig_fn)

   with open(uniquejobid+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")
   return

def send_batchjob_to_local_server(afile,bfile,loopref_checked,pdb_blacklist,uniquejobid):
   outdir = os.path.join(rundir_path,str(uniquejobid))
   looprem_checked = "no";
   lcl_fn = "BatchJob.sh"
   lcl_file_path = os.path.join(rundir_path,uniquejobid,lcl_fn)
   ig_fn = "ignore_list.txt"
   ig_file_path = os.path.join(rundir_path,uniquejobid,ig_fn)

   if pdb_blacklist:
      pdb_blacklist = pdb_blacklist.split(',')
      with open("ignore_list.txt", 'w') as ignore_list_file:
         for ignorepdb in pdb_blacklist:
            ignore_list_file.write(ignorepdb+'\n')

   arecords = list(SeqIO.parse(afile, "fasta"))
   brecords = list(SeqIO.parse(bfile, "fasta"))
         
   with open(lcl_file_path, 'a') as lcl_file:
      lcl_file.write('#!/bin/sh\n')
      prefixnum = 0
      tcrdata['jobid'] = uniquejobid
      tcrdata['tcrmodel'] = []  
      for (arecord,brecord) in zip(arecords, brecords):
         prefixnum += 1
         aregexres = get_vdomain_using_regex(arecord.seq, "A")
         bregexres = get_vdomain_using_regex(brecord.seq, "B")
         tcrdata['tcrmodel'].append({  
               'prefixnum': prefixnum,
               'aseq_user': str(arecord.seq.strip()),
               'bseq_user': str(brecord.seq.strip()),
               'aseq_vdomain': aregexres.groups()[0],
               'bseq_vdomain': bregexres.groups()[0]
         })
         subdir = os.path.join(outdir,str(prefixnum))
         lcl_file.write('mkdir -p '+subdir+'\n')
         lcl_file.write('cd '+subdir+'\n')
         rtcrcommand = "-mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s "  % (arecord.seq,brecord.seq) 
         rtcrcommand += " -out:prefix " + str(prefixnum) + "_ "
         if (loopref_checked == "yes"):
            rtcrcommand += "-refine_tcr_cdr3_loops -loops::max_inner_cycles 50 "
         if pdb_blacklist:
            rtcrcommand += "-ignore_list ignore_list.txt "
         commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/tcrmodel.static.linuxgccrelease " + rtcrcommand +  " > res.out "
         lcl_file.write(commandline+'\n')
         
   with open(uniquejobid+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")

   command = ['bash', lcl_fn]
   p = subprocess.Popen(command)
   return

@app.route('/batchsubmit', methods=['POST', 'GET'])
def batchsubmit():
   if request.method == 'POST':
      if not 'abatchfile' in request.files:
         return render_template("error.html", errormsg="No sequence file uploaded for TCR Alpha chain")
      if not 'bbatchfile' in request.files:
         return render_template("error.html", errormsg="No sequence file uploaded for TCR Beta chain")
      af = request.files['abatchfile']
      bf = request.files['bbatchfile']
      pdb_blacklist = request.form.get('pdbblacklist')
      loopref_checked = "yes" if (request.form.get("lr"))else "no"
      if ( (af and allowed_file(af.filename)) and (bf and allowed_file(bf.filename)) ):
         uniquejobid = create_and_cd_to_unique_dir("TCRM")#TCRM for multiple tcr modeling
         af.save(secure_filename(af.filename))         
         bf.save(secure_filename(bf.filename))         
         arecords = list(SeqIO.parse(af.filename, "fasta"))
         asize = len(arecords)
         brecords = list(SeqIO.parse(bf.filename, "fasta"))
         bsize = len(brecords)
         if not (asize == bsize):
            errormsg = "No. of sequences in the input files are not same"
            return render_template("error.html", errormsg=errormsg)
         #send_batchjob_to_local_server(af.filename,bf.filename,loopref_checked,pdb_blacklist,uniquejobid)
         send_batchjob_to_ibbr_cluster(af.filename,bf.filename,loopref_checked,pdb_blacklist,uniquejobid)
         return redirect(url_for('mtcr', jobid=uniquejobid))
      else:
         return render_template("error.html", errormsg="Error in input file.")


@app.route('/disulfidize_results/<jobid>')
def disulfidize_results(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("disulfidize_error.html", errormsg="Job ID not found")
   os.chdir(outdir)

   resfiles = glob.glob("*_*.pdb")
   if resfiles:
      scorefname = os.path.join(rundir_spath,str(jobid),'score.fasc')
      scorefile = os.path.join(rundir_path,str(jobid),'score.fasc')
      commandline = "grep -v total_score " + scorefile + " | grep -v SEQUENCE | sort -nrk2 | tail -1"
      p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out,err = p.communicate()
      bfile = (out.split(' ')[-1]).strip() + ".pdb"
      modelfname = os.path.join(rundir_spath,str(jobid),bfile)
      return render_template("disulfidize_completed.html", scorefname=scorefname, modelfname=modelfname, resfiles=resfiles, jobid=jobid, bfile=bfile)
   else:
      outfname = os.path.join(rundir_spath,str(jobid),'res.out')
      return render_template("disulfidize_failed.html", outfname=outfname)

'''            
@app.route('/disulfidize')
def disulfidize():
   return render_template("disulfidize.html")
'''

@app.route('/run_disulfidize', methods=['POST', 'GET'])
def run_disulfidize():
   if request.method == 'POST':
      pdbfile = request.files['pdbfile']
      if pdbfile.filename == '':
         return render_template("disulfidize_error.html", errormsg="No file uploaded ")

      uniquedirname = create_and_cd_to_unique_dir("disul")
      #return render_template("error.html", errormsg=uniquedirname)
      pdbfile.save(secure_filename(pdbfile.filename))
      xmlfile = os.path.join(rundir_path,"disulfidize","disulfidize.xml")
      flagfile = os.path.join(rundir_path,"disulfidize","flags")
      commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -database /www/cgi-bin/rosetta/Rosetta/main/database @"+flagfile+" -file:s " + str(pdbfile.filename) +  " -parser:protocol " + str(xmlfile) + " > res.out "
      p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      #Do not call communicate if you do not want to wait for the command to complete
      out,err = p.communicate()
      return redirect(url_for('disulfidize_results', jobid=uniquedirname))
   else:
      return render_template("disulfidize_error.html", errormsg="Error in input file.")

if __name__ == '__main__':
   #app.run(debug = True)
   app.run()
   
