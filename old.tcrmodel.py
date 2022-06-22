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

app = Flask(__name__)

tcrdata = {}


rundir = "tcrmodel_runs"
#Works if rundirpath is one directory down the app (tcrmodel) folder
rundir_path = os.path.abspath(os.path.join(app.root_path, '..', rundir))
#rundir server path
rundir_spath = "/"+rundir

#template db path for pdb and seq
tmplt_db_pdb = "/www/cgi-bin/rosetta/Rosetta/tools/tcr/pdb/"
tmplt_db_seq = "/www/cgi-bin/rosetta/Rosetta/tools/tcr/seq/"

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

def submit_job_on_cluster(uniquejobid,jobscriptile,ig_file_path):
   HOST="cluster.ibbr.umd.edu"
   USER="gowthamanr"
   PASSWD="SaRa1213!"
   remote_dir_path = "/home/gowthamanr/TCRmodeller/" + uniquejobid
   ssh_client = paramiko.SSHClient()
   ssh_client.load_system_host_keys()
   ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
   ssh_client.connect(HOST, username=USER, password=PASSWD)
   sftp = ssh_client.open_sftp()
   sftp.mkdir(remote_dir_path)
   sftp.chdir(remote_dir_path)
   sftp.put(pbs_file_path, jobscriptfile)
   if os.path.isfile(ig_file_path):
      sftp.put(ig_file_path, ig_fn)
   command = "cd %s && sbatch %s" % (remote_dir_path, jobscriptfile)
   ssh_stdin, ssh_stdout, ssh_stderr = ssh_client.exec_command(command)
   return

@app.route('/layout')
def layout():
   return render_template("layout.html")

@app.route('/')
@app.route('/index')
def index():
   return render_template("index.html")

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
   
@app.route('/about')
def about():
   return render_template("about.html")

@app.route('/help')
def help():
   return render_template("help.html")

@app.route('/links')
def links():
   return render_template("links.html")

@app.route('/test')
def test():
   return render_template("test.html")

@app.route('/processjob/<aseq>/<bseq>/<loopref_checked>/<pdb_blacklist>')
def processjob(aseq,bseq,loopref_checked,pdb_blacklist):
   if (aseq and bseq):
    #remove everying but characters from input string
      aseq=''.join(i for i in aseq if i.isalpha())
      bseq=''.join(i for i in bseq if i.isalpha())
      tcrdata['aseq_user']=aseq
      tcrdata['bseq_user']=bseq
   else:
      errormsg = "No input TCR sequence"
      return render_template("error.html", errormsg=errormsg)

   a_regexp_res = get_vdomain_using_regex(aseq,"A") 
   if (not a_regexp_res):
      errormsg = "No TCR Variable domain sequence identified from the Alpha chain input"
      return render_template("error.html", errormsg=errormsg)
   a_vd_seq = a_regexp_res.groups()[0]
   a_gm_seq = a_regexp_res.groups()[1]+a_regexp_res.groups()[2]

   b_regexp_res = get_vdomain_using_regex(bseq,"B") 
   if (not b_regexp_res):
      errormsg = "No TCR variable domain sequence identified from the Beta chain input"
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
      #send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      send_job_to_local_server(uniquedirname,rtcrcommand)
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
   agm_tmplt = seq_match_from_fasta_files_list( a_gm_seq, gm_multidb, pdb_blacklist )
   if agm_tmplt is not None:
      tcrdata['acdr1_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr1_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr1_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      tcrdata['acdr2hv4_tmplt_pdb'] = str(agm_tmplt.id)[:6].upper()
      tcrdata['acdr2hv4_tmplt_pdbid'] = str(agm_tmplt.id)[:4].upper()
      tcrdata['acdr2hv4_tmplt_pdb_chain'] = str(agm_tmplt.id)[5:6].upper()
      agm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(agm_tmplt.id), vd_multidb)
      agm_tmplt_segs = get_cdr_from_seq_by_aho_num(agm_tmplt_vd_seq, "A", True)
      tcrdata['acdr1_tmplt_seq'] = str(agm_tmplt_segs[0])
      tcrdata['acdr2hv4_tmplt_seq'] = str(agm_tmplt_segs[2])
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
      bgm_tmplt_vd_seq = seq_by_id_from_fasta_files_list(str(bgm_tmplt.id), vd_multidb)
      bgm_tmplt_segs = get_cdr_from_seq_by_aho_num(bgm_tmplt_vd_seq, "B", True)
      tcrdata['bcdr1_tmplt_seq'] = str(bgm_tmplt_segs[0])
      tcrdata['bcdr2hv4_tmplt_seq'] = str(bgm_tmplt_segs[2])
   else:
      multidb = glob.glob(tmplt_db_seq+"*_TCR_FW.fasta")
      find_templates(bseq_fw, multidb, "bfw", "Beta Framework", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR1_*.fasta")
      find_templates(bseq_cdr1, multidb, "bcdr1", "Beta CDR1", None, pdb_blacklist)
      multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR2HV4_*.fasta")
      find_templates(bseq_cdr2hv4, multidb, "bcdr2hv4", "Beta CCDr2 & HV4", None, pdb_blacklist)

   #CDR3
   cdr3_multidb = glob.glob(tmplt_db_seq+"*_TCR_CDR3_*.fasta")
   find_templates(aseq_cdr3, cdr3_multidb, "acdr3", "Alpha CDR3", None, pdb_blacklist)
   find_templates(bseq_cdr3, cdr3_multidb, "bcdr3", "Beta CDR3", None, pdb_blacklist)
   #override cdr3 seq with cdr3_extnd seq
   cdr3a_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['acdr3_tmplt_id'], vd_multidb)
   print  "cdr3a_tmplt_vd_seq", cdr3a_tmplt_vd_seq
   cdr3a_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3a_tmplt_vd_seq, None, True)
   tcrdata['acdr3_tmplt_seq'] = str(cdr3a_tmplt_segs[6])
   cdr3b_tmplt_vd_seq = seq_by_id_from_fasta_files_list(tcrdata['bcdr3_tmplt_id'], vd_multidb)
   cdr3b_tmplt_segs = get_cdr_from_seq_by_aho_num(cdr3b_tmplt_vd_seq, None, True)
   tcrdata['bcdr3_tmplt_seq'] = str(cdr3b_tmplt_segs[6])

   with open(uniquedirname+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")
   return redirect(url_for('rtcr', jobid=uniquedirname))

def send_job_to_local_server(uniquedirname,rtcrcommand):
   commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/tcr.static.linuxgccrelease -database /www/cgi-bin/rosetta/Rosetta/main/database " + rtcrcommand +  " > res.out "
   p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   #Do not call communicate if you do not want to wait for the command to complete
   #out,err = p.communicate()
   return 

def send_job_to_ibbr_cluster(uniquejobid,rtcrcommand):
   jobscriptfile = "RosettaTCR.slurm"
   pbs_file_path = os.path.join(rundir_path,uniquejobid,jobscriptfile)
   ig_fn = "ignore_list.txt"
   ig_file_path = os.path.join(rundir_path,uniquejobid,ig_fn)

   with open(pbs_file_path, 'a') as pbs_file:
      pbs_file.write('#!/bin/sh\n')
      pbs_file.write('#SBATCH --export=ALL\n')
      pbs_file.write('#SBATCH -J RosettaTCR\n')
      pbs_file.write('#SBATCH -N 1, -n 16, -t 72:00:00,mem=512M\n')
      #pbs_err_file = uniquejobid + ".err"
      #pbs_out_file = uniquejobid + ".out"
      #pbs_file.write('#PBS -o %s\n' % pbs_out_file)
      #pbs_file.write('#PBS -e %s\n' % pbs_err_file)
      pbs_file.write('cd $SLURM_SUBMIT_DIR\n')
      pbs_file.write('rundir=${PWD##*/}\n')
      commandline = "/home/gowthamanr/Rosetta/main/source/bin/tcr.linuxgccrelease " + rtcrcommand +  " > res.out "
      pbs_file.write(commandline+'\n')
      pbs_file.write("if test -f res.out ; then scp res.out gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/res.out ;fi\n")
      pbs_file.write("chmod 777 $SLURM_SUBMIT_DIR/tcrmodel.pdb\n")
      pbs_file.write("if test -f tcrmodel.pdb ; then scp tcrmodel.pdb gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcrmodel.pdb ;fi\n")
      pbs_file.write("if ! test -f tcrmodel.pdb ; then echo ERROR > tcr.fail; scp tcr.fail gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcr.fail ;fi\n")

   submit_job_on_cluster(uniquejobid,jobscriptfile,ig_file_path)
   return 

@app.route('/submitjob', methods = ['POST', 'GET'])
def submitjob():
   aseq = request.form['alphachain']
   bseq = request.form['betachain']
   if not (aseq or bseq):
      return render_template("error.html", errormsg="No TCR sequence entered!")
   #simil_cutoff = request.form.get('simcutoff')
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   #looprem_checked = "yes" if (request.form.get("lr"))else "no"
   #uniquedirname = processjob(aseq,bseq,simil_cutoff,loopref_checked,looprem_checked,pdb_blacklist)
   return redirect(url_for('processjob', aseq=aseq,bseq=bseq,loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

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
   if os.path.isfile(outpdb):
      renumber_tcrpdb_by_aho_number(outpdb)
      add_ss_header_to_pdbfile(outpdb)
      outpdb_spath = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
      stdout = get_pdb_templates(tcrjsondata)
      return render_template("rtcrcompleted.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   
   elif os.path.isfile(errorfile):
      return render_template("rtcrfailed.html", rtcrjobid = jobid)   
   else:
      return render_template("rtcrqueued.html", jobid=jobid, tj=tcrjsondata)   

@app.route('/mtcr/<jobid>')
def mtcr(jobid):
   outdir = os.path.join(rundir_path,str(jobid))
   if not os.path.exists(outdir):
      return render_template("error.html", errormsg="Job ID not found")
   os.chdir(outdir)
   errorfile = os.path.join(rundir_path,str(jobid),'tcr.fail')
   pdblist = glob.glob("*/*_tcrmodel.pdb")
   if pdblist:
      for pdbfile in pdblist:
         outpdb = os.path.join(rundir_path,str(jobid),pdbfile)
         renumber_tcrpdb_by_aho_number(outpdb)
         add_ss_header_to_pdbfile(outpdb)
      tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
      with open(tcrjsonfile) as json_data:
         tcrjsondata = json.load(json_data)
      #zip files
      zipf = zipfile.ZipFile(str(jobid)+'.zip', 'w', zipfile.ZIP_DEFLATED)
      for pdbfile in pdblist:
         zipf.write(pdbfile)
      zipf.close()
      #return render_template("error.html", errormsg="Job ID not found")
      return render_template("mtcrcompleted.html", mtcrjobid = jobid, tj=tcrjsondata,rundir_spath=rundir_spath)
   elif os.path.isfile(errorfile):
      return render_template("rtcrfailed.html", rtcrjobid = jobid)   
   else:
      return render_template("mtcrqueued.html", jobid=jobid)   

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
   return render_template("rtcrcompleted.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   

@app.route('/viewmodel/<jobid>')
def viewmodel(jobid):
   modelfname = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
   tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
   with open(tcrjsonfile) as json_data:
      tcrjsondata = json.load(json_data)
   return render_template("view.html", jobid=jobid, modelfname=modelfname, tj=tcrjsondata)

@app.route('/mviewmodel/<jobid>/<prefixnum>')
def mviewmodel(jobid,prefixnum):
   modelfname = os.path.join(rundir_spath,str(jobid),str(prefixnum),str(prefixnum)+'_tcrmodel.pdb')
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
   else:
      return render_template("error.html", errormsg="Job ID not found")

@app.route('/testsubmit', methods=['POST', 'GET'])
def testsubmit():
   if request.form.get("aspecies"):
      aspecies = request.form['aspecies']
   else:
      return render_template("error.html", errormsg="Species not selected!")
   if request.form.get("bspecies"):
      bspecies = request.form['bspecies']
   else:
      return render_template("error.html", errormsg="Species not selected!")
   if request.form.get("trav"):
      trav = request.form['trav']
   else:
      return render_template("error.html", errormsg="TRAV gene not selected!")
   if request.form.get("traj"):
      traj = request.form['traj']
   else:
      return render_template("error.html", errormsg="TRAJ gene not selected!")
   if request.form.get("trbv"):
      trbv = request.form['trbv']
   else:
      return render_template("error.html", errormsg="TRBV gene not selected!")
   if request.form.get("trbj"):
      trbj = request.form['trbj']
   else:
      return render_template("error.html", errormsg="TRBJ gene not selected!")
   if request.form.get("acdr"):
      acdr = request.form['acdr']
   else:
      return render_template("error.html", errormsg="CDR3 sequence for alpha chain not entered!")
   if request.form.get("bcdr"):
      bcdr = request.form['bcdr']
   else:
      return render_template("error.html", errormsg="CDR3 sequence for beta chain not entered!")
   pdb_blacklist = request.form.get('pdbblacklist')
   loopref_checked = "yes" if (request.form.get("lr"))else "no"
   aseq =  trav+acdr.upper()+traj
   bseq =  trbv+bcdr.upper()+trbj
   return redirect(url_for('processjob', aseq=aseq,bseq=bseq,loopref_checked=loopref_checked,pdb_blacklist=pdb_blacklist))

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
   jobscriptfile = "RosettaTCR.pbs"
   pbs_file_path = os.path.join(rundir_path,uniquejobid,jobscriptfile)
   ig_fn = "ignore_list.txt"
   ig_file_path = os.path.join(rundir_path,uniquejobid,ig_fn)

   if pdb_blacklist:
      pdb_blacklist = pdb_blacklist.split(',')
      with open("ignore_list.txt", 'w') as ignore_list_file:
         for ignorepdb in pdb_blacklist:
            ignore_list_file.write(ignorepdb+'\n')
            
   arecords = list(SeqIO.parse(afile, "fasta"))
   brecords = list(SeqIO.parse(bfile, "fasta"))
         
   with open(pbs_file_path, 'a') as pbs_file:
      pbs_file.write('#!/bin/sh\n')
      pbs_file.write('#PBS -V\n')
      pbs_file.write('#PBS -S /bin/sh\n')
      pbs_file.write('#PBS -N RosettaTCR\n')
      pbs_file.write('#PBS -l nodes=1:ppn=16,walltime=72:00:00,mem=512mb\n')
    #pbs_file.write('#PBS -q default\n')
      pbs_err_file = uniquejobid + ".err"
      pbs_out_file = uniquejobid + ".out"
      #pbs_file.write('#PBS -o %s\n' % pbs_out_file)
      #pbs_file.write('#PBS -e %s\n' % pbs_err_file)
      pbs_file.write('cd $PBS_O_WORKDIR\n')
      pbs_file.write('rundir=${PWD##*/}\n')
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
         pbs_file.write('mkdir -p $PBS_O_WORKDIR/'+str(prefixnum)+'\n')
         pbs_file.write('cd $PBS_O_WORKDIR/'+str(prefixnum)+'\n')
         rtcrcommand = "-mute all -ignore_zero_occupancy false -renumber_pdb -per_chain_renumbering -alpha %s -beta %s "  % (arecord.seq,brecord.seq) 
         rtcrcommand += " -out:prefix " + str(prefixnum) + "_ "
         if (loopref_checked == "yes"):
            rtcrcommand += "-refine_tcr_cdr3_loops "
         if pdb_blacklist:
            rtcrcommand += "-ignore_list ignore_list.txt "
         commandline = "/home/gowthamanr/Rosetta/main/source/bin/tcr.linuxgccrelease " + rtcrcommand +  " > res.out "
         pbs_file.write(commandline+'\n')
         
      pbs_file.write("chmod -R 777 $PBS_O_WORKDIR\n")
      pbs_file.write("if test -f $PBS_O_WORKDIR/1/1_tcrmodel.pdb ; then scp -r $PBS_O_WORKDIR/*/ gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir ;fi\n")
      pbs_file.write("if ! test -f $PBS_O_WORKDIR/1/1_tcrmodel.pdb ; then echo ERROR > tcr.fail; scp tcr.fail gowthamanr@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcr.fail ;fi\n")

   submit_job_on_cluster(uniquejobid,jobscriptfile,ig_file_path)

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
         commandline = "/www/cgi-bin/rosetta/Rosetta/main/source/bin/tcr.static.linuxgccrelease " + rtcrcommand +  " > res.out "
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
      af = request.files['abatchfile']
      bf = request.files['bbatchfile']
      pdb_blacklist = request.form.get('pdbblacklist')
      loopref_checked = "yes" if (request.form.get("lr"))else "no"

      if af.filename == '':
         return render_template("error.html", errormsg="No sequence file uploaded for TCR Alpha chain")
      if bf.filename == '':
         return render_template("error.html", errormsg="No sequence file uploaded for TCR Beta chain")
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
         send_batchjob_to_local_server(af.filename,bf.filename,loopref_checked,pdb_blacklist,uniquejobid)
         #send_batchjob_to_ibbr_cluster(af.filename,bf.filename,loopref_checked,pdb_blacklist,uniquejobid)
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
            
@app.route('/disulfidize')
def disulfidize():
   return render_template("disulfidize.html")

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
   
