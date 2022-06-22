#!/bin/bash

echo "Weekly PDB updates:"
mv pdb_seqres.txt.gz pdb_seqres.txt.gz.old
mv pdb_seqres.txt pdb_seqres.txt.old
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gunzip -c pdb_seqres.txt.gz > pdb_seqres.txt


echo "Updated PDB sequences and structures"
mv tcr_hmmer/hmmer.hits.pdbseqres.alpha.out tcr_hmmer/hmmer.hits.pdbseqres.alpha.out.old
mv tcr_hmmer/hmmer.hits.pdbseqres.beta.out tcr_hmmer/hmmer.hits.pdbseqres.beta.out.old
hmmsearch tcr_hmmer/tcr.alpha.hmm pdb_seqres.txt > tcr_hmmer/hmmer.pdbseqres.alpha.out
hmmsearch tcr_hmmer/tcr.beta.hmm pdb_seqres.txt > tcr_hmmer/hmmer.pdbseqres.beta.out

perl tcr_hmmer/filter_hmm_output.pl tcr_hmmer/hmmer.pdbseqres.alpha.out 1e-23 > tcr_hmmer/hmmer.hits.pdbseqres.alpha.out
perl tcr_hmmer/filter_hmm_output.pl tcr_hmmer/hmmer.pdbseqres.beta.out 1e-23 > tcr_hmmer/hmmer.hits.pdbseqres.beta.out

awk {'print $9'} tcr_hmmer/hmmer.hits.pdbseqres.alpha.out > tcr_hmmer/TCRalpha_pdbid.txt
awk {'print $9'} tcr_hmmer/hmmer.hits.pdbseqres.beta.out > tcr_hmmer/TCRbeta_pdbid.txt
cat  missing_tcr_beta.txt >> tcr_hmmer/TCRbeta_pdbid.txt

#remove pdbids that already exist in the /www/cgi-bin/rosetta/Rosetta/main/database/additional_protocol_data/tcr/pdb/
bash remove_old_pdbids.sh

mv tcr_alpha.fasta tcr_alpha.fasta.old
mv tcr_beta.fasta tcr_beta.fasta.old

python parse_tcr_domains_from_pdbseq.py 
if [ -s "tcr_alpha.fasta" ] && [ -s "tcr_beta.fasta" ]
then
    acount=$(grep -c "^>" tcr_alpha.fasta)
    bcount=$(grep -c "^>" tcr_beta.fasta)
    echo
    echo "Updated TCR template library: $acount alpha domain and $bcount beta domain structures"
    echo
    echo "Downloading new structures..."
    python get_structures_from_fasta_file.py tcr_alpha.fasta
    python get_structures_from_fasta_file.py tcr_beta.fasta
    echo
    echo "Done!"
    echo
fi
