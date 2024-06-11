import anndata as ad
import numpy as np
import scanpy as sc
import sys
sys.path.append("../scripts")
from hic_functions import hic_zoom, translocation_contact_finder
from filter_scHiC import get_cells, cis_ratio, do_plot
from CUTAG_lib.cutag.utilities.clustering import wanted_leiden
import pandas as pd

outdir="../results/1_exploratory"
anndata_object=sys.argv[1]

#################################################################

#Set parameters

reso = 100_000
trans_reso = 1_000_000
min_count_per_cell = 5000
too_close = 0 # 10_000
min_counts_feature = 7
demux_file="../data/dipc_files/demultiplexed_ditags_fixed_demu.tsv"

#################################################################

#Get good cells

print("Getting cells...")

cells = get_cells(demux_file)
cell_ratio, cell_total, ratio_cut, total_cut, good_cells = cis_ratio( cells, ratio_cut=0.55, total_cut=min_count_per_cell)

#################################################################

#Open bias file

fh = open ("../data/dipc_files/biases_final_dataset.tsv")
biases_cis = {}
biases_trans = {}
for line in fh:
    c, p1, _, v = line.split()
    p1 = int(p1)
    v = float(v)
    biases_cis[c, p1 // reso] = v
    biases_trans.setdefault((c, p1 // trans_reso), []).append(v)
    
for p in biases_trans:
    biases_trans[p] = np.nanmedian(biases_trans[p])

#################################################################

# Initialize an empty dictionary
my_dict = {}

# Open and read the file with two fields
with open('../data/cell_tags.txt', 'r') as file:
    for line in file:
        # Split each line into key and value (assuming space as the separator)
        key, value = line.strip().split()
        
        # Assign the value as the key and the key as the value in the dictionary
        my_dict[value] = key
        
#################################################################

relapse_cells = {key: value for key, value in my_dict.items() if int(value)< 193}
diagnose_cells = {key: value for key, value in my_dict.items() if int(value)>= 193}

translocation_file_pd = pd.read_csv('../data/Miltelman_BALL_translocations.csv', delimiter=',', header=0, names=['BALL',  'translocation', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'refs'])


adts = ad.read_h5ad(anndata_object)

translocation = "t(1;5)(p22;p22)"

print("Getting hic map plot ...")
cells2 = hic_zoom(relapse_cells, good_cells, demux_file, 1000000, 1000000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, outdir)

print("Getting UMAP with translocation...")
translocation_contact_finder(adts, "chr1", 123400000, 125100000, "chr5", 48800001, 59600000, translocation, trans_reso, cells2)

sc.pl.umap(adts, color=['group',translocation], size=200, save=f"_{translocation}_embedding.svg")