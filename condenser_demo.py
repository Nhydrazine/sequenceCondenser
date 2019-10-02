import sys;
import numpy as np, pandas as pd;
import sequenceCondenser as sc;
################################################################################
# This is just a demonstration of the process used by the sequence condenser.  #
# Feel free to modify this script and experiment, that's what this script is   #
# intended for.                                                                #
################################################################################
print("Sequence condenser demonstration.");
print();

print("Original sequence strings:");
seqs = pd.Series({
	'a':'AAAAAAAA',
	'b':'AAAAAAAB',
	'c':'CCCCCAAA',
});
print(seqs);
print();

print("Split into DataFrame:");
seqs = sc.splitSequences(seqs);
print(seqs);
print();

print("Dot-consensus of all sequences:");
print( ''.join( sc.forgiven_consensus(seqs) ) );
print("100% conserved residues are retained, but mismatch positions are dots.");
print();

print("Remove all 100% conserved residues to speed up processing:");
seqs = sc.nonConserved(seqs);
print(seqs);
print();

print("Dot-consensus of all sequences:");
print( ''.join( sc.forgiven_consensus(seqs) ) );
print("Consensus is all dots because we removed conserved residues.");
print();

print("Count non-dot mismatches between first sequence and all others:");
print( sc.mismatches_to_ref( seqs.iloc[0,:], seqs.iloc[1:,:] ) );
print();

print("Identify sequences that are within ft=3 of first sequence:");
print("(within 3 non-dot mismatches to first sequence)");
ft_match = sc.condense_to_ref( seqs.iloc[0,:], seqs.iloc[1:,:], 3 );
print( ft_match.tolist() );
print();

print("Get the dot-consensus of those sequences:");
dot_cons = sc.forgiven_consensus( seqs.loc[ ft_match,: ] );
print(''.join(dot_cons));
print();

print("Replace those sequences with their dot-consensus, call it 'cluster 1'.");
print("So now the sequences are:");
seqs = seqs.drop( ft_match );
dot_cons.name = 'cluster 1';
seqs = seqs.append( dot_cons );
print(seqs);
print();

print("This is the core process that is repeated.");
print("When dot-consensuses are compared, the dot residues don't");
print("count as mismatches.");
print();
print("sequenceCondenser.condense() is the process described above.");
print();
print("sequenceCondenser.condense_multipass() repeats the process");
print("until the dot-consensuses can no longer be clustered together");
print("at the specified forgiveness threhsold.");
print();
print("sequenceCondenser.condense_phase() repeats the multipass");
print("process across increasing forgiveness thresholds until either");
print("all sequences cluster into a single group, or until a maximum");
print("forgiveness threshold is reached.");
print();
print("End of demo.");
