# sequenceCondenser

**Author** Nicholas E. Webb

## Description
The sequence condenser is a python-based sequence clustering algorithm used to explore sequence relationships across a range of relatedness.  The sequence condenser uses a **forgiven consensus**, also referred to as a **dot consensus**, that is used to forgive sequence mismatches when clustering.  The number of forgiven mismatches is defined by the **forgiveness threshold** (`ft`).  With `ft=0` only exact sequence matches will cluster together.  At increasing `ft>0` sequences cluster in progressively larger *nearest-neighbor-like* groups, resembling the physical process of condensation.  At some point, `ft` will be large enough to group all sequences into a single cluster.  The result is a progressive, dynamic view of sequence relationships across an increasing stringency threshold.  This can help emphasize both small and large sequence differences that may be associated with some phenotype.

**Note:** The clusters formed by this algorithm will sometimes differ depending on the order in which the sequences are processed.  This may occur, for example, when a rare variant sequence is the first comparative reference.  For simplicity, the order of sequences in the original FASTA file is preserved by the algorithm, so that running the same file will result in the same clusters.  This property reflects the fact that the algorithm exclusively works with *relative* differences between sequences.

**Note:** This algorithm was specifically designed to explore relationships across discontinuous regions of HIV envelope sequences, particularly the variable loops.  Through this use, we've found that the clusters assigned are sensitive to alignment ambiguity.  For example, the position of a segment that is surrounded by deletions.  The sequence `---RLTNE---` will be clustered separately from `--RLTNE----`.  This algorithm is actually really good a highlighting these trouble spots.

## Quick Start

### Components
* **sequenceCondenser.py**: the sequence condenser module.
* **condenser.py**: runs phase condenser for sequences in a FASTA file (must be aligned!).
* **render.py**: renders cluster/phenotype/sequence display as pdf.
* **condenser_demo.py**: demonstrates the condensation process with text output.
* **example.sh**: bash script that runs *condenser.py* using sample sequences and alignments in *samples* folder.  Also serves as a demonstration of command line arguments for *condenser.py*.
* **samples/sample_sequences.fasta**: small subset of aligned infant/maternal HIV Env quasispecies sequences used by *example.sh*.
* **samples/sample_phenotypes.csv**: sample phenotype file used by *example.sh*.
* **samples/HIV_env-alignment.csv**: sample alignment nomenclature for *samples/sample_sequences.fasta* alignment.

### Requirements

You will need python 3.6+ and the numpy, pandas, argparse, packages.  If you want to use the renderer you will need the matplotlib package (including matplotlib.gridspec).

You will also need a way to execute python scripts with command-line arguments.

You will also need aligned sequences in FASTA format and, optional CSV files containing phenotype values and/or alignment nomenclatures (for the renderer).

### Example Script
You can use **example.sh** and it's breakdown below to start running analyses.  **example.sh** is a sample bash script that demonstrates the usage and arguments of *condenser.py* and *render.py*.  It executes *condenser.py* as follows:
```
python condenser.py \
	samples/sample_sequences.fasta \
	--maxft 40 \
	--out sample_condenser_output.csv
```
providing arguments for the sequence FASTA file, the max ft cutoff and the name of the output file.

Then it runs the *render.py* script demonstrating all of its arguments:
```
python render.py \
	sample_condenser_output.csv \
	samples/sample_sequences.fasta \
	--alignment samples/HIV_env_alignment.csv \
	--region 130,469 \
	--phenotype samples/sample_phenotypes.csv \
	--out sample_renderer_output.pdf \
	--width 11 \
	-nonconserved
```
providing the required condenser output file and sequence file arguments with the following optional arguments:

* **alignment**: a CSV file with one column and no headers.  The column is a sequential list of position *names* according to some alignment nomenclature.  The length of this list must equal the length of the aligned sequences.  These *names* are used in the sequence display axis.  If omitted, positions are referred to using a sequential numeric index starting at zero.

* **region**: a two-element comma separated start and end position for sequence display.  The start and end positions *must be consistent with the nomenclature used* for referring to positions.  If not specified, the entirety of all sequences will be rendered.  This can take a while.

* **phenotype**: an optional CSV file containing three columns: the sequence name, a numeric phenotype value and a numeric error bar value.  Missing values are not plotted.  Only the phenotypes of sequences in the condenser output are used, missing/undefined phenotypes are okay.

* **out**: the name of an output file for the rendered pdf, *must include pdf extension* (default is renderer_output.pdf).

* **width**: the width of the pdf in inches (default 11).

* **nonconserved**: is a switch, when present, only non-conserved residues will be included in the sequence display.  The default is to render all residues (conserved or not), which can be time consuming.  Note that conserved positions are not informative to this analysis.

## The Algorithm

### Single Pass Condenser

* for each sequence in a list of sequences
  * choose the current sequence as a *reference*
  * for all other sequences in the list
    * count the number of mismatches to *reference*
    * *(note that `.` does not count as a mismatch)*
    * if the number of mismatches is <= `ft`
      * assign to cluster with *reference*
      * pop current sequence off list to prevent further comparison
    * else
      * continue
* for each unique cluster assignment
  * generate a **dot consensus** or **forgiven consensus** of all sequences in the cluster by replacing mismatch characters with `.` character
  * *(all the sequences now match and the characteristic pattern of `.` characters can be used to identify the original sequences in the cluster)*
  * reduce or **condense** the original sequence list to a list of **forgiven consensuses** keyed by arbitrary cluster ID
* return **condensed** list of consensuses

### Multi-Pass Condenser
Repeats the single pass algorithm until no more **condensation** occurs (e.g. the number of consensuses coming out is the same as the number going in).  This results in an exhaustive condensation of the sequences.

### Phase Condenser
Repeats the multi-pass algorithm across increasing `ft` until all sequences fall into the same cluster, or until a maximum `ft` threshold is reached.  This results in a view of sequence relationships across a range of mismatch thresholds.

## The Module
The sequenceCondenser module source is fully documented.  For examples, check out condenser_demo.py, which runs through a demo of the single pass condenser.  

## Phase Condenser Script
A complete phase-condenser script is provided (condenser.py).  Check the source for documentation.  This script accepts FASTA files for sequence input, and generates a CSV output where each column is an `ft`, and each row is a sequence.  The value in each cell is the cluster id for the sequence at the indicated `ft`.  The output is sorted to make cluster assignments more obvious.

## Renderer
I've included a python script that renders condenser output with phenotype values and sequence display.  Check the source for documentation and use.
