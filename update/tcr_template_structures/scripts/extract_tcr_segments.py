#!/usr/bin/env python


import sys, os, re
import gzip
from sys import argv
from Bio import SeqIO

from TCR_functions import *

if (len(sys.argv) < 2):
   print ("Usage: python extract_tcr_segments.py <tcrseq> <TAG: EX: A or B for alpha or beta chain>")
   sys.exit()

script, inseq, TAG = argv

#segs = assign_CDRs_using_REGEX(inseq, TAG)
segs = assign_CDRs_using_REGEX_with_cdr3extnd(inseq, TAG)
print segs.groups()
