import sys, os, argparse;
import numpy as np, pandas as pd;
import sequenceCondenser as sc;
"""Runs full condenser algorithm on a FASTA file of sequences.

Uses ``sequenceCondenser.py`` to cluster sequences together across increasing
forgiviess thresholds (``ft``).  The clusters, across all evaluated ``ft``,
along with the original sequences are output into a CSV-formatted file.

Arguments
---------
sequence_file : str
	The name of the FASTA-formatted file containing the sequences to be
	evaluated.

--out : str
	(optional) The name of the CSV file output to.  Defaults to
	``condenser_output.csv``.

--maxft : int
	(optional) The maximum ``ft`` that will abort the run when reached.

Output
------
	This script outputs a CSV-formatted file whose rows are indexed by the
	sequence ID/name (first column), followed by columns that are indexed
	by ``ft`` number.  The values in these columns are the cluster IDs that
	each sequence was assigned to at the given ``ft``.  The last column is
	headed ``sequence``, and contains the corresponding sequence as a single
	string of characters.

Notes
-----
	The sequence file must be in FASTA format.  This was done to avoid
	overloading command-line arguments with additional parameters needed for
	other formats.  The entire FASTA ID line (beginning with ``>``) is taken
	to be the full name of the sequence.

"""
################################################################################
parser = argparse.ArgumentParser();
parser.add_argument(	'sequence_file',
						help =	('File containing sequences.  Can be csv ' +
								 '(.csv), fasta (.fasta), clustal (.clustal) ' +
								 'or sqlite3 compatible database (.sqlite).'),
);
parser.add_argument(	'--out',
						help = 	('Name of output cluster file (csv format). ' +
								 'Default is condenser_output.csv.'),
						default='condenser_output.csv'
);
parser.add_argument(	'--maxft',
						help = ('Maximum cutoff FT value'),
						type=int,
						default=60,
);
args = parser.parse_args();
################################################################################
# LOAD SEQUENCES ###############################################################
################################################################################
base, ext = os.path.splitext(args.sequence_file);
# FASTA only, too much diversity will complicate renderer arguments.
if ext not in ['.fasta']:
	print(	"Invalid file extension for sequence file " +
			str(args.sequence_file) + "...");
	sys.exit();
if not os.path.isfile( args.sequence_file ):
	print("Can't find sequence file " + str(args.sequence_file) + "...");
	sys.exit();
#------------------------------------------------------------------------------#
with open(args.sequence_file,'r') as fh:
	seqs = sc.loadFASTA(fh);
################################################################################
# CONDENSE #####################################################################
################################################################################
# split sequences into pandas.DataFrame, each column is a character.
print("Splitting");
seqs = sc.splitSequences(seqs);

# remove conserved residues, they won't count anyway.
print("Restricting to non-conserved residues");
seqs = sc.nonConserved(seqs);

# run sequence condenser to get a pandas.DataFrame of clusters indexed by
# sequence ID/name, versus FT.
print("Condensing...");
z = sc.condense_phase(seqs, verbose=True, max_ft=args.maxft);

# add sequence strings in the last column.
z['sequence'] = sc.joinSequences(seqs);

# output.
print("Writing to " + args.out + "...");
z.to_csv(args.out);

print("Done.");
