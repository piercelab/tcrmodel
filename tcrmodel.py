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
from TCR_functions import *

app = Flask(__name__)

tcrdata = {}


rundir = "tcrmodel_runs"
#Works if rundirpath is one directory down the app (tcrmodel) folder
rundir_path = os.path.abspath(os.path.join(app.root_path, '..', rundir))
#rundir server path
rundir_spath = "/"+rundir

#template db path for pdb and seq
tmplt_db_pdb = "/home/ragul/Rosetta/main/database/sampling/tcr/pdbs/"
tmplt_db_seq = "/home/ragul/Rosetta/main/database/sampling/tcr/seqs/"

tcrdata['rundir_spath'] = rundir_spath

@app.route('/layout')
def layout():
   return render_template("layout.html")

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
   aregexres = assign_CDRs_using_REGEX_webserver_2version(aseq, "A")
   bregexres = assign_CDRs_using_REGEX_webserver_2version(bseq, "B")

   if (not aregexres):
      errormsg = "No TCR sequence identified from the Alpha chain input"
      return render_template("error.html", errormsg=errormsg)
   if (not bregexres):
      errormsg = "No TCR sequence identified from the Beta chain input"
      return render_template("error.html", errormsg=errormsg)

   #create and cd to unique directory                                                                 
   basename = "TCR"
   random_tag = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(3))
   datetime_tag = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
   uniquedirname = "_".join([("".join([basename, random_tag])), datetime_tag])
   rundir = os.path.join(rundir_path,uniquedirname)


   if not os.path.exists(rundir):
      os.makedirs(rundir)
   os.chmod(rundir, 0777)
   os.chdir(rundir)

   aseqfile = os.path.join(rundir,"aseq.txt")
   bseqfile = os.path.join(rundir,"bseq.txt")
   with open(aseqfile, 'w') as aseq_file:
      for i in range(0, len(aseq), 60):
         aseq_file.write(aseq[i:i+60]+'\n')
   with open(bseqfile, 'w') as bseq_file:
      for i in range(0, len(bseq), 60):
         bseq_file.write(bseq[i:i+60]+'\n')

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
      rtcrcommand += "-refine_tcr_cdr3_loops "
   if ( (loopref_checked == "yes") or (looprem_checked == "yes") ):
      send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      #send_job_to_local_server(uniquedirname,rtcrcommand)
   else:
      send_job_to_ibbr_cluster(uniquedirname,rtcrcommand)
      #send_job_to_local_server(uniquedirname,rtcrcommand)

   #find template for rendering in html
   #create tcr json file with template and sequence details
   tcrdata['jobid'] = uniquedirname
   tcrdata['aseq_vdomain'] = aregexres.groups()[0]
   tcrdata['bseq_vdomain'] = bregexres.groups()[0]
   
   tcrdata['a_fw1'] = aregexres.groups()[1]
   tcrdata['a_fw2'] = aregexres.groups()[3]
   tcrdata['a_fw3'] = aregexres.groups()[5]
   tcrdata['a_fw4'] = aregexres.groups()[7]
   aseq_fw = aregexres.groups()[1] + aregexres.groups()[3] + aregexres.groups()[5] + aregexres.groups()[7]
   tcrdata['aseq_fw'] = aseq_fw
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_FW.fasta")
   afw_tmplt = score_alignment_from_fasta_files_list( aseq_fw, multidb, None, pdb_blacklist )
   if all(afw_tmplt):
      tcrdata['afw_tmplt_pdb'] = str(afw_tmplt[1].id)[:6].upper()
      tcrdata['afw_tmplt_pdbid'] = str(afw_tmplt[1].id)[:4].upper()
      tcrdata['afw_tmplt_seq'] = str(afw_tmplt[1].seq)
   else:
      errormsg = "No Template found for TCR Alpha Framework segment"
      return render_template("error.html", errormsg=errormsg)

   tcrdata['b_fw1'] = bregexres.groups()[1]
   tcrdata['b_fw2'] = bregexres.groups()[3]
   tcrdata['b_fw3'] = bregexres.groups()[5]
   tcrdata['b_fw4'] = bregexres.groups()[7]
   bseq_fw = bregexres.groups()[1] + bregexres.groups()[3] + bregexres.groups()[5] + bregexres.groups()[7]
   tcrdata['bseq_fw'] = bseq_fw
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_FW.fasta")
   bfw_tmplt = score_alignment_from_fasta_files_list( bseq_fw, multidb, None, pdb_blacklist )
   if all(bfw_tmplt):
      tcrdata['bfw_tmplt_pdb'] = str(bfw_tmplt[1].id)[:6].upper()
      tcrdata['bfw_tmplt_pdbid'] = str(bfw_tmplt[1].id)[:4].upper()
      tcrdata['bfw_tmplt_seq'] = str(bfw_tmplt[1].seq)
   else:
      errormsg = "No Template found for TCR Beta Framework segment"
      return render_template("error.html", errormsg=errormsg)

   orientation_template_file = tmplt_db_seq+"TCR1_ORIENTATION.seq"
   ori_tmplt = find_orientation_template(aseq_fw.strip(), bseq_fw.strip(), orientation_template_file)
   if all(ori_tmplt):
      tcrdata['ori_tmplt'] = str(ori_tmplt[0])[:4].upper()
      tcrdata['ori_tmplt_Apdb'] = str(ori_tmplt[0])[:6]
      tcrdata['ori_tmplt_Bpdb'] = str(ori_tmplt[1])[:6]
   else:
      errormsg = "No Template found for orientaion of Alpha and Beta chain sequences"
      return render_template("error.html", errormsg=errormsg)

   aseq_cdr1 = aregexres.groups()[2]
   tcrdata['aseq_cdr1'] = aseq_cdr1
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR1_*.fasta")
   acdr1_tmplt = score_alignment_from_fasta_files_list( aseq_cdr1, multidb, None, pdb_blacklist )
   if all(acdr1_tmplt):
      tcrdata['acdr1_tmplt_pdb'] = str(acdr1_tmplt[1].id)[:6].upper()
      tcrdata['acdr1_tmplt_pdbid'] = str(acdr1_tmplt[1].id)[:4].upper()
      tcrdata['acdr1_tmplt_seq'] = str(acdr1_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR1 of Alpha chain sequence"
      return render_template("error.html", errormsg=errormsg)

   aseq_cdr2hv4 = aregexres.groups()[4]
   tcrdata['aseq_cdr2hv4'] = aseq_cdr2hv4
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR2HV4_*.fasta")
   acdr2hv4_tmplt = score_alignment_from_fasta_files_list( aseq_cdr2hv4, multidb, None, pdb_blacklist )
   if all(acdr2hv4_tmplt):
      tcrdata['acdr2hv4_tmplt_pdb'] = str(acdr2hv4_tmplt[1].id)[:6].upper()
      tcrdata['acdr2hv4_tmplt_pdbid'] = str(acdr2hv4_tmplt[1].id)[:4].upper()
      tcrdata['acdr2hv4_tmplt_seq'] = str(acdr2hv4_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR2HV4 of Alpha chain sequence"
      return render_template("error.html", errormsg=errormsg)

   aseq_cdr3 = aregexres.groups()[6]
   tcrdata['aseq_cdr3'] = aseq_cdr3
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR3_*.fasta")
   acdr3_tmplt = score_alignment_from_fasta_files_list( aseq_cdr3, multidb, None, pdb_blacklist )
   if all(acdr3_tmplt):
      tcrdata['acdr3_tmplt_pdb'] = str(acdr3_tmplt[1].id)[:6].upper()
      tcrdata['acdr3_tmplt_pdbid'] = str(acdr3_tmplt[1].id)[:4].upper()
      tcrdata['acdr3_tmplt_seq'] = str(acdr3_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR3 of Alpha chain sequence"
      return render_template("error.html", errormsg=errormsg)

   bseq_cdr1 = bregexres.groups()[2]
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR1_*.fasta")
   bcdr1_tmplt = score_alignment_from_fasta_files_list( bseq_cdr1, multidb, None, pdb_blacklist )
   tcrdata['bseq_cdr1'] = bseq_cdr1
   if all(bcdr1_tmplt):
      tcrdata['bcdr1_tmplt_pdb'] = str(bcdr1_tmplt[1].id)[:6].upper()
      tcrdata['bcdr1_tmplt_pdbid'] = str(bcdr1_tmplt[1].id)[:4].upper()
      tcrdata['bcdr1_tmplt_seq'] = str(bcdr1_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR1 of Beta chain sequence"
      return render_template("error.html", errormsg=errormsg)

   bseq_cdr2hv4 = bregexres.groups()[4]
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR2HV4_*.fasta")
   bcdr2hv4_tmplt = score_alignment_from_fasta_files_list( bseq_cdr2hv4, multidb, None, pdb_blacklist )
   tcrdata['bseq_cdr2hv4'] = bseq_cdr2hv4
   if all(bcdr2hv4_tmplt):
      tcrdata['bcdr2hv4_tmplt_pdb'] = str(bcdr2hv4_tmplt[1].id)[:6].upper()
      tcrdata['bcdr2hv4_tmplt_pdbid'] = str(bcdr2hv4_tmplt[1].id)[:4].upper()
      tcrdata['bcdr2hv4_tmplt_seq'] = str(bcdr2hv4_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR2HV4 of Beta chain sequence"
      return render_template("error.html", errormsg=errormsg)

   bseq_cdr3 = bregexres.groups()[6]
   multidb = glob.glob(tmplt_db_seq+"*_TCR1_CDR3_*.fasta")
   bcdr3_tmplt = score_alignment_from_fasta_files_list( bseq_cdr3, multidb, None, pdb_blacklist )
   tcrdata['bseq_cdr3'] = bseq_cdr3
   if all(bcdr3_tmplt):
      tcrdata['bcdr3_tmplt_pdb'] = str(bcdr3_tmplt[1].id)[:6].upper()
      tcrdata['bcdr3_tmplt_pdbid'] = str(bcdr3_tmplt[1].id)[:4].upper()
      tcrdata['bcdr3_tmplt_seq'] = str(bcdr3_tmplt[1].seq)
   else:
      errormsg = "No Template found for CDR3 of Beta chain sequence"
      return render_template("error.html", errormsg=errormsg)

   with open(uniquedirname+'.json', 'w') as tcr_json_file:
      js = json.dumps(tcrdata)
      tcr_json_file.write(js+"\n")
   return redirect(url_for('rtcr', jobid=uniquedirname))

def send_job_to_local_server(uniquedirname,rtcrcommand):
   commandline = "/home/ragul/Rosetta/main/source/bin/tcr.linuxgccrelease -database /home/ragul/Rosetta/main/database " + rtcrcommand +  " > res.out "
   p = subprocess.Popen(commandline, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out,err = p.communicate()
   return

def send_job_to_ibbr_cluster(uniquejobid,rtcrcommand):
   pbs_fn = "RosettaTCR.pbs"
   pbs_file_path = os.path.join(rundir_path,uniquejobid,pbs_fn)
   ig_fn = "ignore_list.txt"
   ig_file_path = os.path.join(rundir_path,uniquejobid,ig_fn)

   with open(pbs_file_path, 'a') as pbs_file:
      pbs_file.write('#!/bin/sh\n')
      pbs_file.write('#PBS -V\n')
      pbs_file.write('#PBS -S /bin/sh\n')
      pbs_file.write('#PBS -N RosettaTCR\n')
      pbs_file.write('#PBS -l nodes=1:ppn=1,walltime=72:00:00,mem=512mb\n')
    #pbs_file.write('#PBS -q default\n')
      pbs_err_file = uniquejobid + ".err"
      pbs_out_file = uniquejobid + ".out"
      #pbs_file.write('#PBS -o %s\n' % pbs_out_file)
      #pbs_file.write('#PBS -e %s\n' % pbs_err_file)
      pbs_file.write('cd $PBS_O_WORKDIR\n')
      pbs_file.write('rundir=${PWD##*/}\n')
      commandline = "/home/gowthamanr/Rosetta/main/source/bin/tcr.linuxgccrelease " + rtcrcommand +  " > res.out "
      pbs_file.write(commandline+'\n')
      pbs_file.write("if test -f res.out ; then scp res.out ragul@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/res.out ;fi\n")
      pbs_file.write("if test -f tcrmodel.pdb ; then scp tcrmodel.pdb ragul@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcrmodel.pdb ;fi\n")
      pbs_file.write("if ! test -f tcrmodel.pdb ; then echo ERROR > tcr.fail; scp tcr.fail ragul@piercelab.ibbr.umd.edu:"+rundir_path+"/$rundir/tcr.fail ;fi\n")

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
   sftp.put(pbs_file_path, pbs_fn)
   if os.path.isfile(ig_file_path):
      sftp.put(ig_file_path, ig_fn)
   command = "cd %s && qsub %s" % (remote_dir_path, pbs_fn)
   ssh_stdin, ssh_stdout, ssh_stderr = ssh_client.exec_command(command)
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
   subprocess.call(["mkdssp", outpdb, dsspfile])
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
    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    cmd = subprocess.Popen([profit_program, '-f', infile, '-h', reference, mobile], stdout=subprocess.PIPE)
    return

def get_pdb_templates(tcrjsondata):
   ref = "tcrmodel.pdb"
   refchainA = "A"
   refchainB = "B"
   #cdr1a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr1_tmplt_pdb'][:4].lower()+tcrjsondata['acdr1_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr1_tmplt_pdb'][5:6]
   acdr1_outpdb = tcrjsondata['acdr1_tmplt_pdb'] + "_Acdr1_tmplt.pdb"
   infile = "Acdr1_profit.in"
   create_profit_infile_cdr1(infile,refchainA,mobchain,acdr1_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr2a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr2hv4_tmplt_pdb'][:4].lower()+tcrjsondata['acdr2hv4_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr2hv4_tmplt_pdb'][5:6]
   acdr2hv4_outpdb = tcrjsondata['acdr2hv4_tmplt_pdb'] + "_Acdr2hv4_tmplt.pdb"
   infile = "Acdr2hv4_profit.in"
   create_profit_infile_cdr2hv4(infile,refchainA,mobchain,acdr2hv4_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr3a
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['acdr3_tmplt_pdb'][:4].lower()+tcrjsondata['acdr3_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['acdr3_tmplt_pdb'][5:6]
   acdr3_outpdb = tcrjsondata['acdr3_tmplt_pdb'] + "_Acdr3_tmplt.pdb"
   infile = "Acdr3_profit.in"
   create_profit_infile_cdr3(infile,refchainA,mobchain,acdr3_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr1b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr1_tmplt_pdb'][:4].lower()+tcrjsondata['bcdr1_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr1_tmplt_pdb'][5:6]
   bcdr1_outpdb = tcrjsondata['bcdr1_tmplt_pdb'] + "_Bcdr1_tmplt.pdb"
   infile = "Bcdr1_profit.in"
   create_profit_infile_cdr1(infile,refchainB,mobchain,bcdr1_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr2b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr2hv4_tmplt_pdb'][:4].lower()+tcrjsondata['bcdr2hv4_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr2hv4_tmplt_pdb'][5:6]
   bcdr2hv4_outpdb = tcrjsondata['bcdr2hv4_tmplt_pdb'] + "_Bcdr2hv4_tmplt.pdb"
   infile = "Bcdr2hv4_profit.in"
   create_profit_infile_cdr2hv4(infile,refchainB,mobchain,bcdr2hv4_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #cdr3b
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['bcdr3_tmplt_pdb'][:4].lower()+tcrjsondata['bcdr3_tmplt_pdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['bcdr3_tmplt_pdb'][5:6]
   bcdr3_outpdb = tcrjsondata['bcdr3_tmplt_pdb'] + "_Bcdr3_tmplt.pdb"
   infile = "Bcdr3_profit.in"
   create_profit_infile_cdr3(infile,refchainB,mobchain,bcdr3_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #orientation templates
   #orientation A
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['ori_tmplt_Apdb'][:4].lower()+tcrjsondata['ori_tmplt_Apdb'][4:6]+"_?_aho.pdb.gz")[0]
   mobchain = tcrjsondata['ori_tmplt_Apdb'][5:6]
   oria_outpdb = tcrjsondata['ori_tmplt_Apdb'] + "_oriA_tmplt.pdb"
   infile = "oriA_profit.in"
   create_profit_infile_ori(infile,refchainA,mobchain,oria_outpdb)
   align_tmplt_and_target(infile,ref,mob)
   #orientation B
   mob = glob.glob(tmplt_db_pdb+tcrjsondata['ori_tmplt_Bpdb'][:4].lower()+tcrjsondata['ori_tmplt_Bpdb'][4:6]+"_?_aho.pdb.gz")[0]
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
      get_pdb_templates(tcrjsondata)
      return render_template("rtcrcompleted.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)   
   elif os.path.isfile(errorfile):
      return render_template("rtcrfailed.html", rtcrjobid = jobid)   
   else:
      return render_template("rtcrqueued.html", jobid=jobid, tj=tcrjsondata)   

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


@app.route('/searchid', methods=['POST', 'GET'])
def searchid():
   if not request.form.get("jobidinput"):
      return render_template("error.html", errormsg="Job ID not entered")
   else:
      jobid = request.form['jobidinput']
      outdir = os.path.join(rundir_path,str(jobid))
      if not os.path.exists(outdir):
         return render_template("error.html", errormsg="Job ID not found")
      outpdb = os.path.join(rundir_path,str(jobid),'tcrmodel.pdb')
      errorfile = os.path.join(rundir_path,str(jobid),'tcr.fail')
      tcrjsonfile = os.path.join(rundir_path,str(jobid),str(jobid)+'.json')
      if not os.path.isfile(tcrjsonfile):
         return render_template("error.html", errormsg="Results not found")
      with open(tcrjsonfile) as json_data:
         tcrjsondata = json.load(json_data)
      if os.path.isfile(outpdb):
         outpdb_spath = os.path.join(rundir_spath,str(jobid),'tcrmodel.pdb')
         return render_template("rtcrcompleted.html", rtcrjobid = jobid, modelfname=outpdb_spath, tj=tcrjsondata)
      elif os.path.isfile(errorfile):
         return render_template("rtcrfailed.html", rtcrjobid = jobid)
      else:
         return redirect(url_for('rtcr', jobid=jobid))

      
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

if __name__ == '__main__':
   #app.run(debug = True)
   app.run()
   
