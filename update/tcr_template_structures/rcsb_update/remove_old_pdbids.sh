#!/bin/bash

ls /www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr/pdb/*.gz > processed_pdbid.txt
sed "s/\/www\/cgi-bin\/rosetta\/Rosetta\/main\/database\/additional_protocol_data\/tcr\/pdb\///g" processed_pdbid.txt | sed "s/_aho.pdb.gz//g" > tmp; mv tmp processed_pdbid.txt
cp tcr_hmmer/TCRbeta_pdbid.txt tmp; for i in $(cat processed_pdbid.txt); do grep -v $i tmp > tmp2; mv tmp2 tmp;done; mv tmp tcr_hmmer/TCRbeta_pdbid.txt;
cp tcr_hmmer/TCRalpha_pdbid.txt tmp; for i in $(cat processed_pdbid.txt); do grep -v $i tmp > tmp2; mv tmp2 tmp;done; mv tmp tcr_hmmer/TCRalpha_pdbid.txt;

