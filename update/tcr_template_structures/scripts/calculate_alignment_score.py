#!/usr/bin/env python

import sys, os
from sys import argv

from scoring_matrices import *
from TCR_functions import *
    
script, template_seq, target_seq = argv

scoring_matrix = "BLOSUM62"
score = score_alignment(str(template_seq).upper(), str(target_seq).upper(), scoring_matrix)
print template_seq, target_seq, score
