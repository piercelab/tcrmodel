#!/usr/bin/env python

import os
from Bio import SeqIO
from TCR_functions import *

#fasta_file =  "tcr_alpha.fasta"
fasta_file =  "tcr_beta.fasta"

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for record in fasta_sequences:   
    if check_tcr_alpha_beta_domain(record.seq):
        print record.id, "Found Alpha and beta domain"
        tcrdomainseq = get_tcr_domain_seq(record.seq,"B")
