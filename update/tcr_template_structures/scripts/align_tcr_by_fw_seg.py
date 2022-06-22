#!/usr/bin/env python                                                                                                                                          

import os
from sys import argv
import subprocess
from subprocess import Popen, PIPE

script, refpdb, mobpdb = argv

def align_using_framework_segment(refpdb, mobpdb):
    profit_infile ="profit_framework.in"
    fittedfile = os.path.basename(mobpdb) + ".fitted.pdb"
    if os.path.exists(fittedfile):
        os.remove(fittedfile)
    f = open(profit_infile,"w")
    f.write("ATOMS N,CA,C,O\n")
    f.write("ZONE CLEAR\n")
    f.write("ZONE -23:-23\n")
    f.write("ZONE 43-55:43-55\n")
    f.write("ZONE 91-106:91-106\n")
    f.write("ZONE 139-:139-\n")
    f.write("FIT\n")
    f.write("WRITE " + fittedfile + "\n")

    profit_program = "/TCRmodeller/programs/ProFitV3.1/src/profit"
    myoutput = open("profit_outfile.out", 'w') 
    cmd = subprocess.Popen([profit_program, '-f', profit_infile, '-h', refpdb, mobpdb], stdout=myoutput)

print mobpdb
align_using_framework_segment(refpdb, mobpdb)
