''''

In this script we generate the clustering of the Dip-C cells, outputing the embedding, cells in each cluster and anndata object. We can specify different parameters such as resolution...

'''

#########################################################################################

# IMPORT PACKAGES

import sys
sys.path.append("../scripts/")

import anndata as ad
import networkx as nx
import numpy as np
import os
import scanpy as sc
from collections import defaultdict
from CUTAG_lib.cutag.utilities.clustering import wanted_leiden
from filter_scHiC import get_cells, cis_ratio, do_plot
from matplotlib import pyplot as plt
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import davies_bouldin_score, calinski_harabasz_score, silhouette_score
import warnings
from hic_functions import process_demux_file_multi, sc_contact_finder

# Disable warnings
warnings.filterwarnings("ignore")


#########################################################################################


# Extract command line arguments

min_count_per_cell = int(sys.argv[1])
reso = int(sys.argv[2])
trans_reso = reso *10
too_close = 10_000
min_cells_count = int(sys.argv[3])
mode = sys.argv[4]
demux_file = "../data/dipc_files/demultiplexed_ditags_fixed_demu.tsv"
only_cis = str(sys.argv[5])

# Check mode; diag must go only with only_cis

if mode not in ["raw", "binarize_1", "bias", "diag"]:
    
    print("Please choose a valid normlaization method: raw, binarize_1 or bias")
    sys.exit(1)
    

""" 
bias
raw
log2(raw)
sigmoid transformation lambda x: (1/(1+e^{-x}))
binarize (0 or 1, 1 being at least 3 interactions in the cell)
binarize (0 or 1, 1 being at least 1 interactions in the cell) """

print(f"Cis Resolution: {reso}")
print(f"Trans Resolution: {trans_reso}")
print(f"Too close threshold: {too_close}")
print(f"Min cells that must contain a given feature: {min_cells_count}")
print(f"Normalization mode: {mode}")


cells = get_cells(demux_file)
cell_ratio, cell_total, ratio_cut, total_cut, good_cells = cis_ratio( cells, ratio_cut=0.55, total_cut=min_count_per_cell )

#########################################################################################

with open('../data/cell_tags.txt', 'r') as file:
    cell_class_dict = dict(line.strip().split() for line in file)

cell_class_dict = {v: k for k, v in cell_class_dict.items()}

# Create a new dictionary with keys having values higher than 193
relapse_cells = {key: value for key, value in cell_class_dict.items() if int(value)< 193}
diagnose_cells = {key: value for key, value in cell_class_dict.items() if int(value)>= 193}
    
###############################################################################################

# Create cell object with features

print("Getting counts...")
cells = {}
fh = open(demux_file)
all_features = set()

cchar = 0
for line in fh:
    if line.startswith('#'):
        cchar += len(line)
        continue
    break
fh.seek(cchar)

for line in fh:
    _, c1, b1, s1, l1, _, _, c2, b2, s2, l2, _, _, ct = line.split()
    # caca, crom1, pos1, caca, caca, caca, caca, crom2, pos2, caca, caca, caca, caca, cell-tag
    if ct not in good_cells:
        continue
    if only_cis == "True":
        if c1 != c2:
            continue
    if c1 == "chrY" or c2 == "chrY":
        continue
    b1, l1, b2, l2 = map(int, [b1, l1, b2, l2])
    # define exact position
    p1 = b1 + l1 / 2
    #print(p1)
    #print(p2)
    p2 = b2 + l2 / 2
    # remove diagonal
    if c1 == c2:
        if abs(int(p1) - int(p2)) < too_close:
            continue
        p1 = int(int(p1) // reso)
        p2 = int(int(p2) // reso)
    else:
        p1 = int(int(p1) // trans_reso)
        p2 = int(int(p2) // trans_reso)
    feature = f"{c1}:{p1}_{c2}:{p2}"
    all_features.add(feature)
    try:
        cells[ct][feature] += 1
    except KeyError:
        try:
            cells[ct][feature] = 1
        except KeyError:
            cells[ct] = {feature: 1}

#####################################################################


#Apply biases

if mode=="bias":
    print("Applying biases")
    #Establish biases
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

    for ct, cv in cells.items():
        
        for f, v in cv.items():
            p1, p2 = f.split('_')
            c1, p1 = p1.split(':')
            c2, p2 = p2.split(':')
            if c1 != c2:
                b1 = biases_trans[c1, int(p1)]
                b2 = biases_trans[c2, int(p2)]
            else:
                b1 = biases_cis[c1, int(p1)]
                b2 = biases_cis[c2, int(p2)]
            v = np.log2(v / b1 / b2)
            cv[f] = 0 if np.isnan(v) else v

                
if mode == "binarize_1":
    print("Applying binarization with 1 count")
    for ct, cv in cells.items():
        
        for f, v in cv.items():
            if v >= 1:
                cv[f] = 1
            else:
                
                cv[f] = 0 

if mode == "raw":
    print("Not normalizing data")


#######################################################################

# Generate sparse matrix

feature_idx = dict((f, i) for i, f in enumerate(all_features))
X = sparse.lil_matrix((len(cells), len(all_features)), dtype=np.double)
for i, (c, ff) in enumerate(cells.items()):
    for f, v in ff.items():       
        X[i, feature_idx[f]] = v
        

#######################################################################

# Normalize by diagonal

if mode == "diag":


    X_normalized_diag = sparse.lil_matrix((len(cells), len(all_features)), dtype=np.double)

    for i_cell, (cell, contact_value) in enumerate(cells.items()):
    
        #print(i_cell, cell, contact_value.keys())
    
        # Initialize a defaultdict to accumulate diagonal sums per chromosome
        chromosome_diagonal_sums = defaultdict(lambda: defaultdict(int))

        # Iterate to get all diagonal sums per chromosome
        for feature, value in contact_value.items():
        # Split feature name at underscore to separate bin numbers and chromosome
            chrom1, chrom2 = feature.split('_')[0].split(':')[0], feature.split('_')[1].split(':')[0]
            bin1 = int(feature.split('_')[0].split(':')[1])
            bin2 = int(feature.split('_')[1].split(':')[1])
        
            # Take absolute difference between bin numbers to get diagonal offset
            diagonal_offset = abs(bin2 - bin1)

            # Accumulate the value in the corresponding diagonal and chromosome
            chromosome_diagonal_sums[chrom1][diagonal_offset] += value
            #Divide matrix value by the sum of the correponding diagonal
            if value != X[i_cell, feature_idx[feature]]:
                print("oh,oh....")
        #Iterarte over features to normalize by the sum of the corresponding diagonal
    
        for feature, value in contact_value.items():

            chrom1, chrom2 = feature.split('_')[0].split(':')[0], feature.split('_')[1].split(':')[0]
            bin1 = int(feature.split('_')[0].split(':')[1])
            bin2 = int(feature.split('_')[1].split(':')[1])
        
            # Take absolute difference between bin numbers to get diagonal offset
            diagonal_offset = abs(bin2 - bin1)
            
            X_normalized_diag[i_cell, feature_idx[feature]] = value / (chromosome_diagonal_sums[chrom1][diagonal_offset])    
            
#######################################################################

# Create Anndata

print("Creating AnnData...")


if mode == "diag":
    
    adts = ad.AnnData(X=X_normalized_diag.tocsc(), obs=list(cells.keys()), 

                  var=np.array(list(feature_idx.keys())), 
                  dtype=X.dtype)
else:
    
    adts = ad.AnnData(X=X.tocsc(), obs=list(cells.keys()), 
                  var=np.array(list(feature_idx.keys())), 
                  dtype=X.dtype)


adts.var_names = np.array(list(feature_idx.keys()))
adts.obs_names = list(cells.keys())
del(adts.obs[0])
del(adts.var[0])

########################################################################

########################################################################

del X, fh

for k in ['cis', 'cis 10kb', 'trans', 'cis ratio', 'cis ratio 10kb']:
    adts.obs[k] = [good_cells[c][k] for c in cells]       # Here cells doesn't have the chromosome 1 features, but good cells does

adts.obs['total'] = adts.obs["cis"] + adts.obs["trans"]

# Get unique cell ids

unique_ids = adts.obs.index.tolist()

# Import data for MT counts

import pandas as pd

adts.obs['sample'] = ''

##################################################################################

# Get cell classes

my_dict = {}

# Open and read the file with two fields
with open('../data/cell_tags.txt', 'r') as file:
    for line in file:
        # Split each line into key and value (assuming space as the separator)
        key, value = line.strip().split()
        
        # Assign the value as the key and the key as the value in the dictionary
        my_dict[value] = key

# Loop through the obs_names in the AnnData object
for cell_name in adts.obs_names:
    # Check if the cell_name is in the dictionary
    if cell_name in my_dict:
        # If there's a match, assign the value from the dictionary to the new column
        adts.obs['sample'][cell_name] = my_dict[cell_name]

adts.obs['group'] = ['relapse' if int(value) < 193 else 'diagnose' for value in adts.obs['sample']]


#####################################################################################

# Filter features and regress out variables

print("Filtering and regressing out...")

sc.pp.filter_genes(adts, min_cells=min_cells_count, inplace=True)

if only_cis == "True":
    sc.pp.regress_out(adts, keys=['cis'])
else:
    sc.pp.regress_out(adts, keys=['trans'])
sc.tl.pca(adts, n_comps=20)
sc.pl.pca_variance_ratio(adts, n_pcs=20, save=f"{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis.svg")


######################################################################################

print("Getting PCs correlated with total counts...")

top_correlated_pca = (0, 0)
pcs_to_remove = []
for v in range(0, 19):
    corr = spearmanr(adts.obsm['X_pca'][:,v], adts.obs['cis'])[0]
    if abs(corr) > top_correlated_pca[0]:
        top_correlated_pca = abs(corr), v + 1   #Set correlation and PC number
    print(f"PCA{v+1:<2} corr. with total count: {corr:.2f}")
    if abs(corr) > 0.40:
        print(f'Removing PC {top_correlated_pca[1]}')
        pcs_to_remove.append(v+1)
        
#######################################################################################

variance_explained = adts.uns['pca']['variance_ratio']

adts.obsm['X_pca'] = np.delete(adts.obsm['X_pca'], pcs_to_remove, axis=1)

#######################################################################################

#Perform clustering

sc.pp.neighbors(adts, n_neighbors=6, n_pcs=7)
adts.obsp['connectivities'].data = np.nan_to_num(adts.obsp['connectivities'].data, copy=False)
wanted_leiden(adts, 2)
sc.tl.umap(adts, min_dist=0.5, spread=2)
color_map = {'diagnose': '#5699ffff', 'relapse': '#fa4777ff'}
sc.pl.umap(adts, color=['group'], size=200, palette=color_map)
#sc.pl.umap(adts, color=['leiden', 'group'], size=200)


sc_contact_finder(adts, "chr1_chr5", cells)

df = adts.obs

# Convert the 'feature_present' column to binary (0 if absent, 1 if present)
df['feature_present_binary'] = df['chr1_chr5'].apply(lambda x: 1 if x > 0 else 0)
# Group the DataFrame by the 'group' column
grouped = df.groupby('group')
# Calculate the count of unique observations with the feature present in each group
unique_feature_present_counts = grouped['feature_present_binary'].sum()
# Calculate the total count of observations in each group
total_counts = grouped.size()
# Calculate the percentage of unique observations with the feature present in each group
feature_present_percentages = (unique_feature_present_counts / total_counts) * 100
# Create a DataFrame to store the results
results_df = pd.DataFrame({
    'group': feature_present_percentages.index,
    'cell with feat': unique_feature_present_counts,
    'feat %': feature_present_percentages.values
})

# Print or return the results
print(f'Feature chr1_chr5')
print(results_df)


sc.tl.rank_genes_groups(adts, "leiden", method="t-test")
sc.pl.rank_genes_groups(adts, n_genes=20, save='ranked_genes.png')


ranked_genes = adts.uns['rank_genes_groups']

# Extract the names of the groups (clusters)
groups = ranked_genes['names'].dtype.names

# Print the marker genes for each cluster
for group in groups:
    print(f"Cluster {group}:")
    print(ranked_genes['names'][group])
    print("\n")

demux_file = "../data/dipc_files/dipc_contacts_overlapping_promoter_regions_with_genes.tsv"
print(f"NUmber of cluster1 differential features: {len(ranked_genes['names']['1'])}")
print(f"NUmber of anndata features: {len(adts.var_names)}")

######################################################################################################################################################

print("Saving results...")

if only_cis == "True":
    umap_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis/umaps/"
    anndata_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis/anndata/"
    clusters_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis/cell_cluster_annotation/"
else:
    umap_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all/umaps/"
    anndata_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all/anndata/"
    clusters_path=f"../results/2_clustering_results/clustering_{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all/cell_cluster_annotation/"

if not os.path.exists(umap_path):
        os.makedirs(umap_path)
        
if not os.path.exists(anndata_path):
        os.makedirs(anndata_path)
        
if not os.path.exists(clusters_path):
    os.makedirs(clusters_path)

first_1000_genes_cluster_0 = ranked_genes['names']["0"][:1000]

if only_cis == "True": 
    save_path=f"{umap_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis_umap.png"
    adts_path=f"{anndata_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis_anndata.hdf5"
    cells_path=f"{clusters_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis_cell_clusters.txt"
    _, _ = process_demux_file_multi(cells, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins=set(first_1000_genes_cluster_0), threshold_bins = set(first_1000_genes_cluster_0), output_file=f"{clusters_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_only_cis_marker_genes_contacts_cluster0.txt")
else:
    save_path=f"{umap_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all_umap.png"
    adts_path=f"{anndata_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all_anndata.hdf5"
    cells_path=f"{clusters_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all_cell_clusters.txt"
    _, _ = process_demux_file_multi(cells, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins=set(first_1000_genes_cluster_0), threshold_bins = set(first_1000_genes_cluster_0), output_file=f"{clusters_path}/{reso}_{min_count_per_cell}_{min_cells_count}_{mode}_all_marker_genes_contacts_cluster0.txt")


plt.savefig(save_path)
adts.write(adts_path)
adts.obs[['leiden']].to_csv(cells_path)

file_name=f"{reso}_{min_cells_count}_{min_count_per_cell}_leiden.svg"
sc.pl.umap(adts, color=['leiden'], size=200, save=file_name)

sc.pl.umap(adts, color=['total'], size=200, save=file_name)

adts_copy = adts.copy()  # Make a copy to avoid modifying the original object
sc.pp.neighbors(adts_copy, n_neighbors=6, n_pcs=7)
adts_copy.obsp['connectivities'].data = np.nan_to_num(adts_copy.obsp['connectivities'].data, copy=False)
wanted_leiden(adts_copy, 2)
sc.tl.umap(adts_copy, min_dist=0.5)
connectivities = adts_copy.obsp['connectivities']


#categorical_values = adts_copy.obs['group']

# Create a graph from the connectivities data
G = nx.Graph()
num_nodes = connectivities.shape[0]
G.add_nodes_from(range(num_nodes))

# Add edges based on the threshold (e.g., 0.1)
threshold = 0.3
for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if connectivities[i, j] > threshold:
                G.add_edge(i, j)

# Get assortativity

cell_types = np.array(adts_copy.obs["group"])

# Create a dictionary to associate cell types with nodes
cell_type_dict = {node: cell_type for node, cell_type in zip(G.nodes(), cell_types)}

# Set cell type as a node attribute
nx.set_node_attributes(G, cell_type_dict, 'cell_type')

# Calculate assortativity coefficient
assortativity_coefficient = nx.attribute_assortativity_coefficient(G, 'cell_type')

print(assortativity_coefficient)