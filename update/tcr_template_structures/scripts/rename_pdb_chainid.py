#!/usr/bin/env python

from sys import argv

script, filename, new_chainid = argv

IN = open(filename)

outfilename = 'tmp_'+filename
OUT = open(outfilename, 'w+')

IN = open(filename)
for line in IN.readlines():
    if line[0:4] == 'ATOM':
        #if line[21:22] == current_chainid:
        line  = line[:21] + new_chainid + line[22:]
        OUT.write(line)

IN.close()
OUT.close()
