#!/usr/bin/env python3
import sys
import subprocess
if (len(sys.argv) < 2):
   print ("Usage: python get_core_pep.py <peptide_seq>")
   sys.exit()

def get_core_binding_pepseq(pseq):
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


#print("Argument List:", str(sys.argv))
pseq=sys.argv[1] 
core_pseq = get_core_binding_pepseq(pseq)
print (core_pseq)
