from Bio import SeqIO

fasta_file = "pdb_seqres.txt" # Input fasta file
wanted_TCRalpha_file = "tcr_hmmer/TCRalpha_pdbid.txt" # Input interesting sequence IDs, one per line
result_TCRalpha_file = "tcr_alpha.fasta" # Output fasta file

wanted_TCRbeta_file = "tcr_hmmer/TCRbeta_pdbid.txt" # Input interesting sequence IDs, one per line
result_TCRbeta_file = "tcr_beta.fasta" # Output fasta file


unwanted_list_of_files = ["tcr_like_structures/CD8_list.txt",  "tcr_like_structures/gd_tcr_list.txt",  "tcr_like_structures/IgNAR_list.txt",  "tcr_like_structures/NITR_list.txt", "tcr_like_structures/mmCIF_list.txt"]
unwanted_seqs = set()
for unwanted_file in unwanted_list_of_files:
    with open(unwanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                unwanted_seqs.add(line)


wanted_TCRalpha = set()
with open(wanted_TCRalpha_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted_TCRalpha.add(line)

wanted_TCRbeta = set()
with open(wanted_TCRbeta_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted_TCRbeta.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_TCRalpha_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in unwanted_seqs:
            continue
        if seq.id in wanted_TCRalpha:
            SeqIO.write([seq], f, "fasta")
            
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_TCRbeta_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in unwanted_seqs:
            continue
        if seq.id in wanted_TCRbeta:
            SeqIO.write([seq], f, "fasta")
