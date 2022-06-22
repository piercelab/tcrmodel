#!/usr/bin/env python

import os
import glob
from TCR_functions import *

afilelist = glob.glob("*.A.aho.pdb")
bfilelist = glob.glob("*.B.aho.pdb")
for afile in afilelist:
    apdbid = afile[:4]
    for matchain in glob.glob(apdbid+"*.B.aho.pdb"):
        dist = calc_alphacys_betacys_distance(afile,matchain)
        if (dist < 25):
            print "complex", dist, afile,matchain
            make_complex(afile,matchain)



