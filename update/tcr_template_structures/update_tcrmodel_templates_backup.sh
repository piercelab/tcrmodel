#!/bin/bash
CWD=${PWD} 
echo $CWD

rm -rf $CWD/updated_templates/*
cd ${CWD}/rcsb_update
${CWD}/rcsb_update/update.sh

cd ${CWD}/parsed_template_files
rm -rf *

${CWD}/template_update_scripts/create_template_files.py ${CWD}/rcsb_update/tcr_alpha.fasta A
${CWD}/template_update_scripts/create_template_files.py ${CWD}/rcsb_update/tcr_beta.fasta B
${CWD}/template_update_scripts/create_orientation_templates.py
${CWD}/template_update_scripts/get_unique_orientation_seq.py
${CWD}/template_update_scripts/create_FW_orientation.py
${CWD}/template_update_scripts/create_GM_orientation.py

#update TCRpMHC templates 
cd ${CWD}/tcrpmhc_templates
${CWD}/template_update_scripts/create_tcrpmhc_templates.sh 

#take backup of old templates
cp -r /www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr /www/tcrmodel/update/tcr_template_structures/updated_templates/previous/tcr
#copy new templates to rosetta db
cp -r /www/tcrmodel/update/tcr_template_structures/updated_templates/tcr /www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr
