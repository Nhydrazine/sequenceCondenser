import numpy as np, pandas as pd;
import os, sys;
#------------------------------------------------------------------------------#
"""Sequence condenser alrogithm and associated functions.

The sequence condenser is used to identify differences between a group of
sequences based on their similarity, through a recursive series of clustering
processes.

The Process
-----------
The sequence condenser clusters sequences together using a forgiveness threshold
(abbreviated ``ft``), which is similar, but not the same, as a hamming distance
or mismatch count.  Once sequences are clustered together into a group, the
entire group is removed and replaced with a dot-consensus sequence.  The
dot-consensus contains all of the sequence residues that are totally conserved
among the cluster, and a period or dot (``.``) at any position that is not
totally conserved among the cluster.  The dot-consensus is then used to
represent all sequences in the cluster, thus, reducing or condensing the
sequences into a smaller number of representative dot-consensus sequences.

The dot-consensus sequences are then clustered in the same way that the
original sequences were, except that the dots do not count as mismatches. The
result is a progressive condensation of sequences into dot-consensuses, and
then a condensation of dot-consensuses.  This repeats until the total number
of dot-consensuses no longer decreases - e.g. no more condensation occurs.
So far, the entire process described occurs under a single forgiveness
threshold, which defines the number of non-dot mismatches allowed among
sequences in a cluster.

Finally, the process described above can be run across increasing forgiveness
thresholds, to observe how the sequences cluster together as one relaxes
the strictness of cluster definitions.

Notes
-----
A sequence Series is a pandas.Series (list-like, or array-like) of
sequential characters in an alphanumeric sequence.

A sequence DataFrame is a pandas.DataFrame of sequence Series'.  The
row indices are the sequence IDs or names, and the column indices can
be used for custom positional nomenclature or alignment nomenclature
on the sequence(s).

A dot-consensus is a sequence that is used to represent a group of other
aligned sequences.  The dot-consensus represents positions that are totally
conserved by the residue that is conserved, and represents positions that are
anything but totally conserved as a period character ``.`` (or dot).

"""
#------------------------------------------------------------------------------#
################################################################################
# CORE SEQUENCE CONDENSER FUNCTIONS ############################################
################################################################################
def mismatches_to_ref(ref,seq):
	#--------------------------------------------------------------------------#
	"""Counts the number of mismatches between one sequence and each individual
	sequence in a group of sequences.

	Counts the number of character mismatches between a single reference
	sequence Series and each individual sequence in a DataFrame of sequence
	Series'. Returns a pandas.Series indexed with the sequencs from the
	DataFrame.  The values in the Series are the number of mismatches between
	the indexed sequence and the reference sequence.  Note that the length of
	the reference Series must be the same as the number of columns in the
	DataFrame.

	Parameters
	----------
	ref : pandas.Series
		A sequence series where each character in the sequence is an
		individual element of the series.
	seq : pandas.DataFrame
		A DataFrame of sequences indexed by ID/name, with each column
		representing a single position of the sequence, e.g. a DataFrame
		of sequence series'.

	Returns
	-------
	pandas.Series
		A series indexed by the sequence IDs/names in ``seq`` whose values
		indicate the number of mismatch positions between the corresponding
		sequence and the ``ref`` sequence.

	"""
	#--------------------------------------------------------------------------#
	# collect all dot positions into a Series and remove those positions from
	# both reference and sequences.
	seq_rmp = seq[seq=='.'].any();
	ref_rmp = ref[ref=='.'];
	rmp = seq_rmp[seq_rmp==True].index.tolist() + ref_rmp.index.tolist();
	seq = seq.drop( columns=rmp );
	ref = ref.drop( rmp );
	# return the mismatches between dot-removed sequences
	return (seq!=ref).sum(axis=1);
################################################################################
def forgiven_consensus(seq):
	#--------------------------------------------------------------------------#
	"""Generates a dot-consensus of sequences in a sequence DataFrame.

	Generate the dot-consneus sequence, as a sequence Series, from a group
	of sequences in a sequence DataFrame.

	Parameters
	----------
	seq : pandas.DataFrame
		A DataFrame of sequences indexed by ID/name, with each column
		representing a single position of the sequence, e.g. a DataFrame
		of sequence series'.

	Returns
	-------
	pandas.Series
		The dot-consensus sequence Series.

	"""
	#--------------------------------------------------------------------------#
	rc = (seq==seq.iloc[0]).sum(axis=0);
	cons = seq.iloc[0].copy();
	cons[ rc[rc<len(seq)].index ] = '.';
	return cons;
################################################################################
def condense_to_ref(ref,seq,ft):
	#--------------------------------------------------------------------------#
	"""Generates a list of all sequences in a sequence DataFrame that are within
	a certain forgiveness threshold (``ft``) of a reference sequence Series.

	Generates a list of all sequences in a sequence DataFrame that are within
	``ft`` mismatches to a reference sequence Series.  The mismatch counts do
	not include potential dot positions when dot-consensus sequences are
	present.

	Parameters
	----------
	ref : pd.Series
		A reference sequence Series to compare to.
	seq : pd.DataFrame
		A sequence DataFrame of sequences.
	ft : int
		The desired forigveness threshold (hamming distance).

	Returns
	-------
	pandas.Index
		A list-like index of sequence IDs/names from the sequence DataFrame
		that are within a forgiveness threshold of ``ft`` to the reference
		sequence Series.

	"""
	#--------------------------------------------------------------------------#
	# return a list of the sequence IDs that match within ft
	mismatch_counts = mismatches_to_ref(ref, seq);	# count non-dot mismatches
	mismatch_counts[ref.name] = 0;
	return mismatch_counts[mismatch_counts<=ft].index;
################################################################################
def condense(seq, ft, max_cycles=100):
	"""Single-pass condenser function that assigns sequences to cluster IDs
	using a forgiveness threshold.

	Assigns sequences to clusters, where each cluster contains sequences with
	``ft`` non-dot mismatches between them.  This is done in a single pass that
	prevents any sequence that is assigned to a cluster, from participating in
	further comparisons.  The resulting clusters may be within ``ft`` non-dot
 	mismatches of one another.

	Parameters
	----------
	seq : pandas.DataFrame
		A sequence DataFrame of sequences.
	ft : int
		The desired forigveness threshold (non-dot mismatches allowed).
	max_cycles : int
		Maximum number of repetition attempts before aborting.

	Returns
	-------
	pandas.Series
		A pandas.Series indexed by sequence ID/name (as in ``seqs`` DataFrame)
		whose values are the cluster IDs

	"""
	#--------------------------------------------------------------------------#
	# build list of all sequence IDs/names and initialize
	seq_queue = seq.index;
	cluster_assignments = pd.Series(0,index=seq.index);
	STOP = 0;		# trigger
	passno = 0;		# current cycle
	# loop until all sequences have been compared (no more in queue)
	while STOP == 0:
		# if no sequences or only one cluster then stop
		if len(seq)<=1: STOP += 1;
		# otherwise condense
		if STOP == 0 and len(seq_queue)>1:
			# get sequences that have not been clustered yet
			lseq = seq.loc[seq_queue];
			# get sequences within ft using first sequence as ref
			cluster = condense_to_ref(lseq.iloc[0], lseq.iloc[1:], ft);
			# assign sequences to cluster ID's by cycle number
			cluster_assignments[cluster] = passno + 1;
			# remove clustered sequences from queue so they aren't compared
			# again
			seq_queue = seq_queue[~seq_queue.isin(cluster)];
		else: STOP += 1;
		if passno > max_cycles:
			raise RuntimeError("Exceeded maximum cycles...");
		passno += 1;
	return cluster_assignments;
################################################################################
def condense_multipass(seq,ft,verbose=False,progress=False,max_cycles=100):
	#--------------------------------------------------------------------------#
	"""Multi-pass condenser that repeatedly clusters dot-consensus sequences
	until no additional clustering occurs at the desgnated forgiveness
	threshold.

	Multi-pass condenser that repeatedly executes the single-pass condenser
	until no more condensation (e.g. clustering) occurs.  The cluster
	assignments returned by the single-pass condenser are converted to dot-
	consensus sequences which are then clustered during the next repetition.
	The first iteration, then, uses the single-pass condenser on the original
	sequences.  All subsequent iterations run the single-pass condenser on the
	dot-consensus sequences of the clusters formed during the previous run.
	Because the dots (``.``) in the dot-consensus are forgiven (do not count
	towards mismatches), the dot-consensus sequences will cluster together in
	ways that the original sequences could not. This condensation, or reduction
	of clusters, continues until no more condensation occurs at the given
	``ft``.

	Parameters
	----------
	seq : pandas.DataFrame
		A sequence DataFrame of sequences.
	ft : int
		The desired forigveness threshold (non-dot mismatches allowed).
	verbose : boolean
	progress : boolean
	max_cycles : int
		Maximum number of repetition attempts before aborting for
		``condense()``.

	Returns
	-------
	pandas.DataFrame
		A sequence pandas.DataFrame indexed by cluster ID, where each sequence
		is the dot-consensus of the cluster.

	"""
	#--------------------------------------------------------------------------#
	useq = seq.copy();
	for z in range(0,40):
		# get number of clusters we're starting with this round
		initial_clusters = len(useq);
		# generate clusters from useq sequences
		clusters = condense(useq, ft, max_cycles=max_cycles );
		# temporary list of cluster consensus sequences
		tseq = [];
		# for every unique cluster id
		for cid in clusters.unique():
			# get list of sequence ids in cluster
			crec = clusters[clusters==cid];
			# generate dot-consensus of the sequences
			c = forgiven_consensus(useq.loc[crec.index]);
			# name the consensus by cluster id
			c.name = cid+1;
			# add the cluster consensus to temporary list
			tseq.append(c);
		if progress: print('.', end="", flush=True);
		# set useq equal to the new consensus sequences and repeat
		useq = pd.DataFrame(tseq, columns=useq.columns);
		# if we got the same number of clusters as we started with
		# no more condensing/clustering can be done
		if initial_clusters == len(useq): break;
		# useq is now a pd.DataFrame of consensus sequences indexed by cluster
		# ID
	return useq;
################################################################################
def condense_phase(seq,verbose=False,max_cycles=100,max_ft=40):
	#--------------------------------------------------------------------------#
	"""Condenses sequences across increasing forgiveness thresholds (``ft``)
	until all sequences fall into the same cluster.

	Progerssively condenses sequences across increasing forgiveness thresholds
	(``ft``) using the multi-pass condenser function at each ``ft``.  This
	process is repeated until either max_ft is reched or until it reaches an
	``ft`` where all sequences fall into the same cluster (maximum forgiveness).
	The maximum forgiveness threhsold depends on the diversity of the sequences
	themselves.

	Parameters
	----------
	seq : pandas.DataFrame
		A sequence DataFrame of sequences.
	verbose : boolean
		Optional.
	max_cycles : int
		Maximum number of repetition attempts before aborting for
		``condense()``.

	Returns
	-------
	pandas.DataFrame
		A DataFrame indexed by sequence ID/name.  Columns are ``ft`` values
		and each cell contains the cluster ID of the sequence at that ``ft``.

	Example
	-------
	>>> z = condense_phase(seqs,progress=True);
	>>> print(z);

	"""
	#--------------------------------------------------------------------------#
	clusters = pd.DataFrame(index=seq.index);
	# index=env, col=ft, val=cluster ID @ ft for all FTs
	for ft in range(0,max_ft+1):
		if verbose==True:
			print(str('[FT ' + str(ft) + '] ').ljust(9), end='', flush=True);
		# get consensus sequences at current FT
		ft_consensuses = condense_multipass(
			seq,
			ft,
			verbose=verbose,
			max_cycles=max_cycles,
		);
		# if we have more than one cluster
		if len(ft_consensuses)>1:
			# one cluster is the terminal condition
			# now get the individual sequence names that fall into each cluster
			# by comparison
			for i,cseq in ft_consensuses.iterrows():
				# trim cluster consensus to only 100% conserved residues
				cons_trim = cseq[cseq!='.'];
				# get all sequences but trim positions where the cluster has '.'
				seqs_trimmed = seq[ cons_trim.index ];
				# now count mismatches between each sequnce and conserved
				# cluster consensus positions
				seqs_match = seqs_trimmed.eq( cons_trim ).sum(axis=1);
				# restrict to only sequences that matched 100% to conserved
				# cluster consnesus (not counting '.')
				seqs_match = seqs_match[seqs_match==len(cons_trim)];
				clusters.loc[seqs_match.index.tolist(),ft] = i;
				if verbose: print('.', end="", flush=True);
		# if only one cluster
		else:
			if verbose: print('.', end="", flush=True);
			clusters.loc[ :,ft ] = 1;
			break;
		print();
	return(clusters);
# MWE:
#z = condense_phase(seqs,progress=True);
#print(z);

################################################################################
# BASIC MANIPULATION FUNCTIONS #################################################
################################################################################
def splitSequences(seqs):
	"""Split a pandas.Series of sequences (``seqs``) into a pandas.DataFrame
	indexed by sequence ID/name, where each sequence position is a single
	column."""
	return pd.DataFrame([ list(r) for r in seqs ], index=seqs.index);

################################################################################
def joinSequences(seqs):
	"""Joins columns in a pandas.DataFrame of sequences (``seqs``) together into
	a pandas.Series of sequence strings."""
	ss = pd.Series([ ''.join(r) for i,r in seqs.iterrows() ]);
	ss.index = seqs.index;
	return ss;

################################################################################
def nonConserved(seqs):
	"""Removes all conserved columns from a sequence pandas.DataFrame."""
	rc = (seqs==seqs.iloc[0]).sum(axis=0);
	return seqs.loc[ :,rc[rc<len(seqs)].index ];

################################################################################
# USEFUL IO/FORMATTING FUNCTIONS ###############################################
################################################################################
def loadFASTA(fh):
	"""Returns pd.Series of sequences as string indexed by sequence id/name
	from file handle pointing to FASTA format sequence file."""
	seqs = {};
	current = "";
	for l in fh.read().split('\n'):
		if l[0]=='>': current = l[1:];
		else:
			if current in seqs.keys():
				seqs[current] += str(l);
			else:
				seqs[current] = str(l);
	return pd.Series(seqs, dtype=str).drop(index=[''], errors='ignore');

def loadCSV(fn,idcol,seqcol):
	"""Loads sequences from a CSV file where idcol is the name of the
	sequence and seqcol is the column of the sequences as strings."""
	return pd.read_csv(	fn,
						header=0,
						index_col=idcol)[seqcol].astype(str);
