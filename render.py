import os, sys, argparse;
import numpy as np, pandas as pd;
import matplotlib.pyplot as plt, matplotlib as mpl;
import matplotlib.gridspec as gridspec;
import sequenceCondenser as sc;
################################################################################
root = os.path.abspath( os.path.dirname( sys.argv[0] ) );
################################################################################
parser = argparse.ArgumentParser();
parser.add_argument(    'condenser_file',
                        help=("sequenceCondenser CSV output file"),
                        type=str,
);
parser.add_argument(    'sequence_file',
                        help=("path/name of FASTA file containing sequences"),
                        type=str,
);
parser.add_argument(    '--phenotypes',
                        help=("(optional) path/name of CSV phenotype file"),
                        default="",
                        type=str,
);
parser.add_argument(    '--alignment',
                        help=("(optional) path/name of CSV file containing "+
                            "sequence alignment nomenclature"),
                        default="",
                        type=str,
);
parser.add_argument(    '--region',
                        help=("(optional) comma separated start,end sequence "+
                            "alignment positions to display (must match "+
                            "nomemclature if specified)"),
                        default="",
                        type=str,
);
parser.add_argument(    '--out',
                        help=("(optional) name of pdf output file, default is "+
                        "renderer_output.pdf"),
                        default="renderer_output.pdf",
);
parser.add_argument(    '--width',
                        help=("(optional) page width in inches, default is 11"),
                        type=float,
                        default=11,
);
parser.add_argument(    '-nonconserved',
                        help=("trim conserved residues from sequence display"),
                        default=False,
                        action='store_true',
);
args = parser.parse_args();
#------------------------------------------------------------------------------#
if not os.path.isfile( args.condenser_file ):
    print("Can't find "+str(args.condenser_file));
    sys.exit();
if not os.path.isfile( args.sequence_file ):
    print("Can't find "+str(args.condenser_file));
    sys.exit();
if args.phenotypes!='' and not os.path.isfile(args.phenotypes):
    print("Can't find "+str(args.phenotypes));
    sys.exit();
if args.alignment!='' and not os.path.isfile(args.alignment):
    print("Can't find "+str(args.alignment));
    sys.exit();
seq_region = args.region.split(',');
if len(seq_region)>2:
    print("Bad definition for --region");
    sys.exit();
#------------------------------------------------------------------------------#
# DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# input_file = os.path.join(root, "condenser_output.csv");
# seq_file = os.path.join(root, "samples", "sample_sequences_baby.fasta");
# align_file = os.path.join(root, "samples", "HIV_env_alignment.csv");
# pheno_file = os.path.join(root, "samples", "sample_phenotypes_baby.csv");
# seq_region = ['138_02','400'];
# #seq_region=[];
# seq_nonconserved = False;
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
input_file = args.condenser_file;
seq_file = args.sequence_file;
pheno_file = args.phenotypes;
align_file = args.alignment;
seq_nonconserved = args.nonconserved;
output_file = args.out;
################################################################################
# LOAD PHENOTYPES
################################################################################
pheno = pd.DataFrame(columns=['sequence','phenotype','error']);
if pheno_file!="":
    if os.path.isfile(pheno_file):
        pheno = pd.read_csv( pheno_file, index_col=0, header=0 );
################################################################################
# LOAD SEQUENCES
################################################################################
base, ext = os.path.splitext(seq_file);
# FASTA only, too much diversity will complicate renderer arguments.
if ext not in ['.fasta']:
	print(	"Invalid file extension for sequence file " +
			str(seq_file) + "...");
	sys.exit();
if not os.path.isfile( seq_file ):
	print("Can't find sequence file " + str(seq_file) + "...");
	sys.exit();
#------------------------------------------------------------------------------#
with open(seq_file,'r') as fh:
	seqs = sc.loadFASTA(fh);
#------------------------------------------------------------------------------#
seqs = sc.splitSequences(seqs);
# not all nomemclature is purely numeric
seqs.columns=seqs.columns.astype(str); # cast str to avoid type conversions
################################################################################
# LOAD ALIGNMENT
################################################################################
align = pd.Series(seqs.columns);
if align_file!=None:
    if os.path.isfile(align_file):
        align = pd.read_csv( align_file, index_col=None, header=None )[0];
try:
    seqs.columns = align;
except Exception as e:
    print("An error occurred while setting alignment:");
    print(str(e));
    sys.exit();
################################################################################
# TRIM SEQUENCE
# (after loading alignment)
################################################################################
# to region if specified
if len(seq_region)==2:
    seqs = seqs.loc[:,seq_region[0]:seq_region[1]];
# trim sequences to non-conserved
if seq_nonconserved:
    seqs = sc.nonConserved(seqs);
################################################################################
# LOAD AND PREPARE CONDENSER OUTPUT
################################################################################
clusters = pd.read_csv( input_file, index_col=0 ).drop(columns=['sequence']);
clusters.columns = clusters.columns.astype(int);
# get informative ft values as series for plot indexing
ftvals = pd.Series( clusters.T.drop_duplicates(keep='last').T.columns );
# subset clusters
clusters = clusters[ftvals];
# sort for visual clarity and joining groups as bars
# sort low to high entropy = reverse ft order
clusters = clusters.sort_values(by=ftvals[::-1].tolist());
# build y values as sequential numerical index
clusters['ix'] = np.arange(0,len(clusters));
# transparent access to sequence name
clusters['name'] = clusters.index.tolist();
################################################################################
# VALIDATION AND FILLING
################################################################################
# check cluster<-->sequence indices (FATAL)
if len(clusters.index[clusters.index.isin(seqs.index)])!=len(clusters):
    missing = clusters.index[~clusters.index.isin(seqs.index)];
    print("Missing sequences for the following: "+
        '.'.join([ str(s) for s in missing ])
    );
    sys.exit();
# check cluster<-->phenotype indices (FILL np.nan)
nopheno = pd.DataFrame(
    columns=pheno.columns,
    index=clusters.index[~clusters.index.isin(pheno.index)]
);
nopheno.loc[nopheno.index,'phenotype'] = np.nan;
nopheno.loc[nopheno.index,'error'] = np.nan;
pheno = pd.concat([pheno, nopheno], ignore_index=False);
################################################################################
# build parameters for each bar
# for every ft: for every unique cluster id:
# ft = ft, group = group, lo = min ix, max = max ix
bars = [];
for ft in ftvals:
    for g in clusters[ft].unique():
        ss = clusters[clusters[ft]==g]; # subset
        bars.append( pd.Series({
            'ft'    : int(ft),          # ft value
            'g'     : g,                # cluster id / group name
            'lo'    : ss['ix'].min(),   # lowest sequence y value
            'hi'    : ss['ix'].max(),   # highest sequence y value
            'ftix'  : int(ftvals[ftvals==ft].index[0]), # ft index number for x
        }));
bars = pd.DataFrame(bars);
################################################################################
# plot page class to help standardize matplotlib graph specs
class plotPage(object):
    def __init__(s,
        width=11,
        height=8.5,
        columns=1,
        rows=1,
        padding_left=0.5,   # applied to each cell
        padding_right=0.3,  # applied to each cell
        padding_top=0.3,    # applied to each cell
        padding_bottom=0.4, # applied to each cell
    ):
        s.width = width;
        s.height = height;
        s.columns = columns;
        s.rows = rows;
        s.cells = [];
        for r in range(0,s.rows):
            s.cells.append([]);
            for c in range(0,s.columns):
                s.cells[r].append({
                    'row'           : r,
                    'column'        : c,
                    'padding_top'   : padding_top,
                    'padding_left'  : padding_left,
                    'padding_right' : padding_right,
                    'padding_bottom': padding_bottom,
                    'width'         : s.width/columns,
                    'height'        : s.height/rows,
                    'left'          : (s.width/columns)*c,
                    'top'           : (s.height/rows)*r,
                });
        s.elements = {
            'header': {
                'padding_top'   : 0.1,
                'padding_left'  : 0,
                'padding_right' : 0,
                'padding_bottom': 0,
                'width'         : s.width,
                'height'        : 0.5,
                'left'          : 0,
                'top'           : 0,
            },
            'footer': {
                'padding_top'   : 0,
                'padding_left'  : 0,
                'padding_right' : 0,
                'padding_bottom': 0.0,
                'width'         : s.width,
                'height'        : 0.5,
                'left'          : 0,
                'top'           : s.height-0.5,
            },
        };
        return;

    def _get_positions(s, r):
        return (
            (r['padding_left']+r['left']) / s.width,
            1-((r['top']+r['height']-r['padding_bottom']) / s.height),
            (r['width']-r['padding_right']-r['padding_left']) / s.width,
            (r['height']-r['padding_bottom']-r['padding_top']) / s.height,
        );

    def get_positions(s, row=0, column=0):
        r = s.cells[row][column];
        return s._get_positions(r);

    def element_positions(s, name):
        return s._get_positions(s.elements[name]);

    def set(s, row=0, column=0, **kwargs):
        props = list(s.cells[0][0].keys());
        for kw in kwargs.keys():
            if kw in props:
                s.cells[row][column][kw] = kwargs[kw];
################################################################################
aa_colors = {
	'D'	: '#ec8677',
	'E'	: '#ec8677',
	'C'	: '#fdec82',
	'M'	: '#fdec82',
	'K'	: '#89b0d0',
	'R'	: '#89b0d0',
	'S'	: '#f3b670',
	'T'	: '#f3b670',
	'F'	: '#90a0c7',
	'Y'	: '#90a0c7',
	'N'	: '#d2eac8',
	'Q'	: '#d2eac8',
	'G'	: '#a7a7a7',
	'L'	: '#bcdc79',
	'V'	: '#bcdc79',
	'I'	: '#bcdc79',
	'A'	: '#808080',
	'W'	: '#f5cfe4',
	'H'	: '#b384b9',
	'P'	: '#f2b06e',
	'X'	: '#555555',
	'-'	: '#ffffff',
	'.'	: '#ff00ff',
	'0'	: '#ff00ff',
};
################################################################################
# determine widths
w_total = args.width;   # total page width in inches
w_per_ft = 0.16;        # ~1/8" ft bars
w_per_res = 0.110;      # sequence residues
w_pheno = 1.6;          # pheno width
w_pad_left = 0.5;   # page pad left
w_pad_right = 0.5;  # page pad right
w_spacing = 0.1;    # space between plots in row
# calculate width of cluster plot
w_clusterplot_yaxlabel = 1;  # padding for y axis labels
w_clusterplot = w_per_ft*len(ftvals)+w_clusterplot_yaxlabel;
# sequence gets remainder
w_sequence = w_total-w_pad_left-w_pad_right-w_pheno-2*w_spacing-w_clusterplot;
w_residues = int(np.floor(w_sequence / w_per_res)); # num residues we can fit
# actual width total:
#print(
#   w_pad_left+w_clusterplot+w_spacing+w_pheno+w_spacing+w_sequence+w_pad_right
#);
################################################################################
# determine heights
h_total = 8;        # landscape
rows = int(np.floor(len(seqs.columns) / w_residues))+1;
h_per_seq = 0.31;   # height per sequence
h_top_spacing = 0.1;    # top padding
h_bot_spacing = 0.5;    # bottom padding includes axis labels
h_row = h_top_spacing + h_per_seq*len(seqs) * h_bot_spacing;
################################################################################
# FIX TODO
# WARNING: plt.subplots returns a 1D array when nrows=1 or ncols=1 and a 2D
# array otherwise.  This causes indexing errors when only one row is needed.
# for now, just make an extra row that will be empty sequence...
if rows<=1: rows=2;
# TODO
# pagination later.  For now, just calculate a total height for all rows
h_total = h_row*rows + 1; # plus some extra padding
page = plotPage(
    width       = w_total,
    height      = h_total,
    columns     = 3,
    rows        = rows,
);
mpl.rcParams['figure.figsize'] = (page.width, page.height);
mpl.rcParams['font.size'] = 8;
mpl.rcParams['font.family'] = 'arial';
# initialize all axes then adjust and link to row/col
fig, ax = plt.subplots(nrows=rows, ncols=page.columns, sharey='row');
cmap = mpl.cm.get_cmap('nipy_spectral');

print("Rendering "+str(rows)+" rows:     ",end='',flush=True);
def show_percent(v):
    p = str(int(np.round(v))).ljust(3)+'%';
    print('\b\b\b\b'+p,end='',flush=True);

for r in range(0,rows):
    show_percent((r/rows)*100);
    page.set(
        column          = 0,
        row             = r,
        left            = w_pad_left,
        top             = r*h_row,
        width           = w_clusterplot,
        height          = h_row,
        padding_left    = w_clusterplot_yaxlabel,
        padding_bottom  = h_bot_spacing,
        padding_top     = h_top_spacing,
        padding_right   = 0,
    );
    page.set(
        column          = 1,
        row             = r,
        left            = w_pad_left+w_clusterplot,
        top             = r*h_row,
        width           = w_pheno+w_spacing,
        height          = h_row,
        padding_left    = w_spacing,
        padding_bottom  = h_bot_spacing,
        padding_top     = h_top_spacing,
        padding_right   = 0,
    );
    page.set(
        column          = 2,
        row             = r,
        left            = w_pad_left+w_clusterplot+w_spacing+w_pheno,
        top             = r*h_row,
        width           = w_sequence,
        height          = h_row,
        padding_left    = w_spacing,
        padding_bottom  = h_bot_spacing,
        padding_top     = h_top_spacing,
        padding_right   = 0,
    );
    for c in range(0,3):
        # set positions
        ax[r,c].set_position( page.get_positions( row=r, column=c ) );
        # specific to first column
        if c==0:
            ax[r,c].set_ylim(-0.5, (len(clusters)-0.5));
            ax[r,c].set_yticks(range(0,len(clusters)));
            ax[r,c].set_yticklabels(clusters.index);
        # specific to rest of columns
        else:
            plt.setp(ax[r,c].get_yticklabels(), visible=False);
        # individual columns
        ax[r,0].set_xlim( -0.5, len(ftvals)-0.5 );
        ax[r,0].set_xticks(range(0,len(ftvals)));
        ax[r,0].set_xticklabels(ftvals.astype(int), rotation=0);
        # column 1: phenotype
        # set limits and labels here if you like
        ax[r,1].grid(
            which='major',
            axis='y',
            color='#000000',
            alpha=0.6,
            linewidth=0.5,
            linestyle='dotted',
        );
        ax[r,2].set_xlim( -0.5, w_residues-0.5);
        ax[r,2].set_xticks( range(0,w_residues) );
        ax[r,2].set_xticklabels(
            seqs.columns[ (w_residues*r):(w_residues*(r+1)) ],
            rotation=90
        );
    # plot bars
    for i,b in bars.iterrows():
        mid = (((b['hi']-b['lo'])/2) + b['lo'] + 0.5) / len(clusters);
        if mid < 0.45:
            textcolor   = '#ffffff';
            boxcolor    = '#ffffff';
        else:
            textcolor   = '#000000';
            boxcolor    = '#ffffff';
        ax[r,0].bar(
            b['ftix'],
            (b['hi']-b['lo'])+1,
            bottom=b['lo']-0.5,
            color=cmap(mid),
            edgecolor=boxcolor,
            linewidth=1,
            width=1,
        );
        ax[r,0].text(
            b['ftix'],
            b['hi']+0.5,
            str(int(b['g'])),
            verticalalignment='top',
            horizontalalignment='center',
            weight='heavy',
            color=textcolor
        )
    # phenotype
    if len(pheno)>0:
        ax[r,1].errorbar(
            pheno.loc[clusters.index,'phenotype'],
            clusters['ix'],
            xerr=pheno.loc[clusters.index,'error'],
            fmt='none',
            elinewidth=1,
            ecolor='#000000',
            capsize=2,
            alpha=0.3
        );
        ax[r,1].scatter(
            pheno.loc[clusters.index,'phenotype'],
            clusters['ix'],
            marker='o',
            color='#888888',
            s=30,
        );
    for t in ax[r,1].get_xticklabels(): t.set_rotation(90);
    # plot residues
    for y in np.arange(0,len(clusters)): # for each sequence
        for x in np.arange(0,w_residues):
            if ((r*w_residues)+x)>=len(seqs.columns): break;
            residue = seqs.loc[
                clusters.index[y],
                seqs.columns[(r*w_residues)+x]
            ];
            try: color = aa_colors[residue];
            except: color='#999999';
            ax[r,2].text(
                x,
                y,
                residue,
                verticalalignment='center',
                horizontalalignment='center',
                bbox=dict(
                    facecolor=color,
                    edgecolor='#000000',
                    pad=1.6,
                    linewidth=0.3
                ),
                fontdict=dict(
                    family='monospace',
                    size=7,
                ),
            );
################################################################################
show_percent(100);
print();
print("Writing to "+str(os.path.join(root,output_file))+"...");
plt.savefig( os.path.join(root,output_file) );
################################################################################
# fin.
