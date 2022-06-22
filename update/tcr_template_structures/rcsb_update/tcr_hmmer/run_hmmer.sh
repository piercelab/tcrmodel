mv hmmer.hits.pdbseqres.alpha.out hmmer.hits.pdbseqres.alpha.out.old
mv hmmer.hits.pdbseqres.beta.out hmmer.hits.pdbseqres.beta.out.old
hmmsearch mmult.alpha.hmm ~/databases/pdb_seqres_9_2015.fa > hmmer.pdbseqres.alpha.out
hmmsearch mmult.beta.hmm ~/databases/pdb_seqres_9_2015.fa > hmmer.pdbseqres.beta.out
hmmsearch ../sequences/trdv.human_mouse.imgt.hmm ~/databases/pdb_seqres_9_2015.fa > hmmer.pdbseqres.delta.out
../scripts/filter_hmm_output.pl hmmer.pdbseqres.delta.out 1e-30 > hmmer.hits.pdbseqres.delta.out
../scripts/filter_hmm_output.pl hmmer.pdbseqres.alpha.out 1e-23 hmmer.hits.pdbseqres.delta.out > hmmer.hits.pdbseqres.alpha.out
../scripts/filter_hmm_output.pl hmmer.pdbseqres.beta.out 1e-23 > hmmer.hits.pdbseqres.beta.out

#diff hmmer.hits.pdbseqres.alpha.out hmmer.hits.pdbseqres.alpha.out.old
#diff hmmer.hits.pdbseqres.beta.out hmmer.hits.pdbseqres.beta.out.old
