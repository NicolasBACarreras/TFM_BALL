'''

This scripts integrates all the lichcic datasets, first with the bootstrap analysis (gets the common contacts to a file), generates the PCA of the pseudobulks and al cell types

'''

#Load modules

import sys
sys.path.append("../scripts/")
import os
import numpy as np
from scipy import sparse
import anndata as ad
from hic_functions import process_demux_file_multi, get_adts_object, jaccard_similarity, process_lichic_file, spearman_correlation_lichic_neo, dict_to_vector, jaccard_similarity_boot, spearman_correlation_multi_lichic
from CUTAG_lib.cutag.utilities.clustering import wanted_leiden
from filter_scHiC import get_cells, cis_ratio, do_plot
from collections import Counter
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import spearmanr
import scanpy as sc
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd

#######################################################################################

# Load anndata object from clustering step

use_adata=False
use_demux=True

if use_adata:
    
    print("Loading anndata...")
    adata = ad.read_h5ad("/home/bsc/bsc059153/PROJECTS/BALL/TFM/results/2_clustering_results/clustering_100000_4000_10_raw_all/anndata/100000_4000_10_raw_all_anndata.hdf5")
    
    class_of_interest = 'diagnose'  # Replace with the class you are interested in
    diagnose_adata = adata[adata.obs['group'] == class_of_interest]
    diagnose_features = diagnose_adata[:, diagnose_adata.X.sum(axis=0) > 0].var_names

    class_of_interest = 'relapse'  # Replace with the class you are interested in
    relapse_adata = adata[adata.obs['group'] == class_of_interest]
    relapse_features = diagnose_adata[:, relapse_adata.X.sum(axis=0) > 0].var_names


# Load demux file
  
if use_demux:
    
    print("Using dip-c file...")
    
    demux_file=sys.argv[1]
    min_count_cell=int(sys.argv[2])
    #min_count_cell=5000
    reso = int(sys.argv[3])  
    #reso = 100000  
    too_close = int(sys.argv[4])     
    #too_close = 10000    

    outdir=f"../results/3_lichic_integration/{min_count_cell}_{reso}_{too_close}_outdir/"
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #Get groups
    
    with open("/home/bsc/bsc059153/PROJECTS/BALL/TFM/data/cell_tags.txt", 'r') as file:

        cell_class_dict = dict(line.strip().split() for line in file)
        cell_class_dict = {v: k for k, v in cell_class_dict.items()}


    relapse_cells = {key: value for key, value in cell_class_dict.items() if int(value)< 193}
    diagnose_cells = {key: value for key, value in cell_class_dict.items() if int(value)>= 193}
    
    #Get good cells
    
    cells = get_cells(demux_file)
    cell_ratio, cell_total, ratio_cut, total_cut, good_cells = cis_ratio(cells, ratio_cut=0.55, total_cut=min_count_cell)
    
    demux_file="/home/bsc/bsc059153/PROJECTS/BALL/TFM/data/dipc_files/dipc_contacts_overlapping_promoter_regions_with_genes.tsv"
    
    cells_relapse, relapse_features = process_demux_file_multi(relapse_cells, good_cells, reso, too_close, demux_file)
    cells_diagnose, diagnose_features = process_demux_file_multi(diagnose_cells, good_cells, reso, too_close, demux_file)
    
    

#######################################################################################

# Load liCHi-C data
lichic_file="/home/bsc/bsc059153/PROJECTS/BALL/TFM/data/lichic_files/BALL_relapse_2_cutoff_5.csv"
lichic_relapse = process_lichic_file(lichic_file, too_close, reso)
red_lines = []
#Filter lichic contacts 

#print(lichic_relapse.keys())
#lichic_relapse = {k: v for (k, v) in lichic_relapse.items() if int(v) > 2 }

#######################################################################################

# Similarity of relapse cell

n_feat_sc, n_feat_lichic, n_feat_inter, similarity_score, common_features_relapse = jaccard_similarity(list(relapse_features), list(lichic_relapse.keys()))

correlation, pvalue = spearman_correlation_lichic_neo(cells_relapse, lichic_relapse)

#Print results to stdout

print(f"There are {n_feat_inter} unique contacts from the lichic data ({n_feat_lichic} total unique contacts) that are present in the schic data ({n_feat_sc} different interactions using relapse data)")
print(f"Similarity score: {similarity_score}")
print(f"Spearman correlation: {correlation}")

red_lines.append(similarity_score) 

# Similarity of diagnose cell

n_feat_sc, n_feat_lichic, n_feat_inter, similarity_score, common_features_diagnose = jaccard_similarity(list(diagnose_features), list(lichic_relapse.keys()))

correlation, pvalue = spearman_correlation_lichic_neo(cells_diagnose, lichic_relapse)

#Print results to stdout

print(f"There are {n_feat_inter} unique contacts from the lichic data ({n_feat_lichic} total unique contacts) that are present in the schic data ({n_feat_sc} different interactions using relapse data)")
print(f"Similarity score: {similarity_score}")
print(f"Spearman correlation: {correlation}")

red_lines.append(similarity_score)

#######################################################################################

# Generate venn diagram

# Get intersection between relapse and diagnose

common_values = set(common_features_diagnose).intersection(common_features_relapse)

# Number of common interactions

print(f"Number of common interactions between sets: {len(list(common_values))}")

# Similarity score BETWEEN sets

similarity_among_sets = len(list(common_values)) / ( len(common_features_diagnose) + len(common_features_relapse) - 2*len(list(common_values)) )

print(f"Similarity between relapse and diagnose: {similarity_among_sets}")

plt.figure(figsize=(10, 6), dpi=100)

# Create a Venn diagram

venn_labels = {'100': 'Common Bins with lichic in Relapse', '010': 'Common Bins in with lichic in Diagnose', '110': 'Intersection'}
venn2(subsets=[set(common_features_relapse), set(common_features_diagnose)], set_labels=('Common Bins in Relapse', 'Common Bins in Diagnose'))

plt.savefig(f"{outdir}venn_diagram.png")

plt.clf()


######################################################################################

# Save contacts

#Note, maybe not just present or not present but maybe proportions

relapse_exclusive_contacts = []
for contact in common_features_relapse:
    if contact not in common_values:
        relapse_exclusive_contacts.append(contact)
        

diagnose_exclusive_contacts = []
for contact in common_features_diagnose:
    if contact not in common_values:
        diagnose_exclusive_contacts.append(contact)
        


#######################################################################################

# Generate output files

print("Saving unique contacts into files...")

_, _, relapse_feat_in_cell = dict_to_vector(cells_relapse)

_, _, diagnose_feat_in_cell = dict_to_vector(cells_diagnose)


count_bins_relapse = {k: v for k, v in relapse_feat_in_cell.items() if v >= 3}
threshold_bins_relapse = set(count_bins_relapse.keys())

count_bins_diagnose = {k: v for k, v in diagnose_feat_in_cell.items() if v >= 3}
threshold_bins_diagnose = set(count_bins_diagnose.keys())

cells_relapse, relapse_features = process_demux_file_multi(relapse_cells, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins = common_features_relapse, threshold_bins = threshold_bins_relapse, output_file="../6_gsea/relapse_exclusive_features.txt")
cells_diagnose, diagnose_features = process_demux_file_multi(diagnose_cells, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins = common_features_diagnose, threshold_bins = threshold_bins_diagnose, output_file="../6_gsea/diagnose_exclusive_features.txt")

""" 
use also all features that overlap with lichic as control

background in enricher chicago promoters or union of single cell relapse and diagnose
"""
 
#######################################################################################

""" # Bootstrapping

print("Performing bootstrapping of both classes...")

# Get counts for diagnose and relapse

    
diagnose_cells_boot, all_diagnose_features = process_demux_file_multi(diagnose_cells, good_cells, reso, too_close, demux_file)
    
relapse_cells_boot, all_relapse_features = process_demux_file_multi(relapse_cells, good_cells, reso, too_close, demux_file)

adts_diagnose = get_adts_object(diagnose_cells_boot, all_diagnose_features)
adts_relapse = get_adts_object(relapse_cells_boot, all_relapse_features)

# Bootstrap diagnose

# Placeholder for results
similarity_scores_diagnose = []

# Number of iterations (bootstrap)
num_iterations = 1000
# Get the number of columns in the AnnData object
columns_in_diagnose = adts_diagnose.shape[1]

for _ in range(num_iterations):

    # Randomly sample n column indices with replacement
    sampled_diagnose_features = np.random.choice(columns_in_diagnose, size=adts_diagnose.n_vars, replace=True)

    # Extract the sampled columns from the AnnData object
    bootstrap_diagnose_adts = adts_diagnose[:, sampled_diagnose_features]

    # Call the similarity score function
    similarity_score = jaccard_similarity_boot(list(bootstrap_diagnose_adts.var_names), list(lichic_relapse.keys()))

    
    # Append the similarity score to the list
    similarity_scores_diagnose.append(similarity_score)

# Print the average similarity score over all iterations
average_similarity_score = np.mean(similarity_scores_diagnose)
print(f"Average similarity score for diagnose bootstrapping: {average_similarity_score}")

#-----------------------------------------------------------------------------------------------------------------

#Bootstrap for relapse

# Placeholder for results
similarity_scores_relapse = []

# Number of iterations (bootstrap)
num_iterations = 1000
# Get the number of columns in the AnnData object
columns_in_relapse = adts_relapse.shape[1]

for _ in range(num_iterations):

    # Randomly sample n column indices with replacement
    sampled_relapse_features = np.random.choice(columns_in_relapse, size=adts_relapse.n_vars, replace=True)

    # Extract the sampled columns from the AnnData object
    bootstrap_relapse_adts = adts_relapse[:, sampled_relapse_features]

    # Call the similarity score function
    similarity_score = jaccard_similarity_boot(list(bootstrap_relapse_adts.var_names), list(lichic_relapse.keys()))

    # Append the similarity score to the list
    similarity_scores_relapse.append(similarity_score)

# Print the average similarity score over all iterations
average_similarity_score = np.mean(similarity_scores_relapse)
print(f"Average similarity score for relapse bootstrapping: {average_similarity_score}")


################


# Plot bootstrapping results

data_range = (min(min(similarity_scores_diagnose), min(similarity_scores_relapse)), max(max(similarity_scores_diagnose), max(similarity_scores_relapse)))
num_bins = 20  # You can adjust the number of bins as needed

plt.figure(figsize=(10, 6), dpi=100)

colors = sns.color_palette("Set2")


#plt.hist(similarity_scores_diagnose, range=data_range, bins=100, color='blue', edgecolor='None', alpha=0.6, label="Diagnose")
#plt.hist(similarity_scores_relapse, range=data_range, bins=100, color='pink', edgecolor='None', alpha=0.6, label="Relapse")

# Add kernel density estimate (KDE) plots as smooth curves

sns.histplot(data=similarity_scores_diagnose, kde=True, color=colors[0], linestyle='None', label="Diagnose")
sns.histplot(data=similarity_scores_relapse, kde=True, color=colors[1], linestyle='None', label="Relapse")
#Add red lines

plt.axvline(x=red_lines[0], color=colors[1], linestyle='--', linewidth=2, label=f'Original relapse similarity value')
plt.axvline(x=red_lines[1], color=colors[0], linestyle='--', linewidth=2, label=f'Original diagnose similarity value')
        

# Set labels and title
plt.xlabel('Similarity scores')
plt.ylabel('Count')
plt.title('Similarity score with bootstrap of diagnose and relapse')

#plt.legend(title="Overlap between significant LIChIC interactions:", fontsize='small')

plt.legend(title="Overlap between significant \n LIChIC interactions:", fontsize='small', bbox_to_anchor=(1,1), loc="upper left", frameon=False)
plt.tight_layout() # para que cuadre una leyenda tan ancha
# Save the plot
plt.savefig(f'{outdir}/bootstrapping_diagnose_relapse.svg', format='svg') """




##########################################################################################################


# Multi lihcic PCA

#########

# Get files
print("Getting lichic files")

directory = '../data/lichic_files/'

# Get all files in the directory
files = os.listdir(directory)

# Filter files that end with 'ibed'
ibed_files = [file for file in files if file.endswith('ibed')]

result_dict = {}

for filename in ibed_files:
    # Extracting the prefix before the first underscore
    prefix = filename.split('_')[0]
    if prefix in ["PC", "immtransB", "Ery", "CLP", "relapse1", "diagnose1"]:
        continue
    full_name =  f"../data/lichic_files/{filename}"
    print(full_name)
    result_dict[prefix] = process_lichic_file(full_name, too_close, reso)

result_dict["relapse2"] = lichic_relapse

#########

unique_contacts = {}
unique_contacts_length = {}
# Iterate over major keys
for cell in result_dict:
    if cell in ["PC", "immtransB","CMP", "Ery"]:
        continue
    # Get contacts for the current cell and store them in a set
    contacts = set(result_dict[cell].keys())
    # Find contacts unique to the current cell
    unique_contacts[cell] = contacts - set().union(*(set(result_dict[key].keys()) for key in result_dict if key != cell))
    unique_contacts_length[cell] = len( contacts - set().union(*(set(result_dict[key].keys()) for key in result_dict if key != cell))) / len(contacts)

total_contacts_per_cell_type = {}
for key, contacts in result_dict.items():
    for contact, value in contacts.items():
        try:
            total_contacts_per_cell_type[key] += value
        except KeyError:
            total_contacts_per_cell_type[key] = value
            
            
#########

plt.bar(unique_contacts_length.keys(), unique_contacts_length.values(), color='skyblue')
plt.xlabel('Cell')
plt.ylabel('Number of Unique Contacts')
plt.title('Number of Unique Contacts per Cell Type')
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Save the plot as a file (e.g., PNG format)
plt.savefig(f'{outdir}/unique_contacts_per_cell_type.png', bbox_inches='tight')

plt.clf()


plt.bar(total_contacts_per_cell_type.keys(), total_contacts_per_cell_type.values(), color='skyblue')
plt.xlabel('Cell')
plt.ylabel('Number of Reads')
plt.title('Number of Reads Contacts per Cell Type')
plt.xticks(rotation=45, ha='right')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Save the plot as a file (e.g., PNG format)
plt.savefig(f'{outdir}/total_contacts_per_cell_type.png', bbox_inches='tight')

plt.clf()


#########


print("Getting anndata object...")

adata = ad.read_h5ad("/home/bsc/bsc059153/PROJECTS/BALL/TFM/results/2_clustering_results/clustering_100000_5000_7_raw_all/anndata/100000_5000_7_raw_all_anndata.hdf5")  #Should be the same as the the parameteres used for the previous comparison. If the adata file doesnt exist run the clustering file first

relapse_cell_tags = []
diagnose_cell_tags = []

for i in adata.obs_names:
    if adata.obs["group"][i] == "relapse":
        relapse_cell_tags.append(i)
    else:
        diagnose_cell_tags.append(i)
        
cells_relapse_pseudo, _ =  process_demux_file_multi(relapse_cells, good_cells, reso, too_close, demux_file)
        
relapse_schic = {}
for key, value in cells_relapse_pseudo.items():
    for contact, values in value.items():
        try:
            relapse_schic[contact] += values
        except KeyError:
            relapse_schic[contact] = values
            
new_relapse_schic = {key: value for key, value in relapse_schic.items() if value > 5}

cells_diagose_pseudo, _ =  process_demux_file_multi(diagnose_cells, good_cells, reso, too_close, demux_file)
        
diagnose_schic = {}
for key, value in cells_diagose_pseudo.items():
    for contact, values in value.items():
        try:
            diagnose_schic[contact] += values
        except KeyError:
            diagnose_schic[contact] = values
            
new_diagnose_schic = {key: value for key, value in diagnose_schic.items() if value > 5}

result_dict["sc_relapse"] = new_relapse_schic
result_dict["sc_diagnose"] = new_diagnose_schic

#########

# Binarize data

for key in result_dict:
    for subkey in result_dict[key]:
        # Assign 1 if value > 0, else leave it at 0
        if result_dict[key][subkey] > 1:
            result_dict[key][subkey] = 1
        else:
            result_dict[key][subkey] = 0

#########
# Get PCA plot

print("Performing PCA of cell types")

# Extract unique keys from all inner dictionaries
inner_keys = sorted(set().union(*[d.keys() for d in result_dict.values()]))

# Create a list of lists to hold the values
feature_idx = dict((f, i) for i, f in enumerate(inner_keys))
X = sparse.lil_matrix((len(result_dict.keys()), len(inner_keys)), dtype=np.double)
for i, (c, ff) in enumerate(result_dict.items()):
    for f, v in ff.items():
        #print(v)       
        X[i, feature_idx[f]] = int(v)


# Convert dictionary of dictionaries to Anndata object
adata_total = ad.AnnData(X=X.tocsc(), obs=list(result_dict.keys()), 
                  var=np.array(list(inner_keys)), 
                  dtype=X.dtype)

adata_total.var_names  = np.array(list(inner_keys))
adata_total.obs_names = list(result_dict.keys())

del(adata_total.obs[0])
del(adata_total.var[0])

sc.pp.filter_genes(adata_total, min_cells=3, inplace=True)

adata_total.obs['total'] = ""
for cell_type in result_dict.keys():    
    adata_total.obs['total'][cell_type] = sum(result_dict[cell_type].values()) 
    
adata_total.obs['total'] = adata_total.obs['total'].astype(float)

sc.pp.regress_out(adata_total, keys=['total'])
sc.tl.pca(adata_total, n_comps=10)

######

#Hierearchrial


adata_total.obsm['X_pca'] = adata_total.obsm['X_pca']  # Ensure PCA is in the AnnData object
linkage_matrix = linkage(adata_total.obsm['X_pca'], method='ward')

# Create a dendrogram plot
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=adata_total.obs_names)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Sample Index')
plt.ylabel('Distance')

# Rotate the x-axis labels to be diagonal
plt.xticks(rotation=45, ha='right')
plt.tight_layout()


plt.savefig(f'{outdir}/hierarchical_clustering_dendrogram.png')
plt.show()

######


labels = adata_total.obs.index  # replace 'label' with the actual column name in your data
plt.figure(figsize=(6, 6))

color_palette = {
    'Mon': 'yellow',
    'HSC': 'black',
    'CMP': 'grey',
    'CLP': 'grey',
    'nCD4': 'lightblue',
    'nCD8': 'lightgreen',
    'nB': 'cyan',
    'memB': 'blue',
    'tB': 'mediumblue',
    'PC': 'lime',
    'GCB': 'darkgreen',
    'immtransB': 'aquamarine',
    'PreB': 'orange',
    'PreProB': 'indianred',
    'ProB': 'firebrick',
    'sc_relapse': 'tomato',
    'sc_diagnose': 'pink',
    'relapse1': 'red',
    'relapse2': 'green',
    'diagnose1': 'purple',
    'Ery': 'olive'
}

# Create a list of colors based on cell types
colors = [color_palette[cell_type] for cell_type in adata_total.obs_names]


# PC1 VS PC2

plt.scatter(adata_total.obsm['X_pca'][:, 0], adata_total.obsm['X_pca'][:, 1], c=colors, s=100)

# Annotate legend with labels
legend_handles = []
for label, color in color_palette.items():
    if label in [ "PC", "immtransB", "CLP", "Ery"]:
        continue
    legend_handles.append(plt.scatter([], [], label=label, color=color))

plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.title('PCA of liCHiC samples')
plt.xlabel(f"PC1 Explained Variability: {adata_total.uns['pca']['variance_ratio'][0]:.3f}")
plt.ylabel(f"PC2 Explained Variability: {adata_total.uns['pca']['variance_ratio'][1]:.3f}")
plt.savefig(f'{outdir}/pca_lichic_cell_types_cmp_preb_gcb_pc1_pc2.svg', format='svg', bbox_inches='tight')
plt.clf()

# PC1 VS PC3

plt.scatter(adata_total.obsm['X_pca'][:, 0], adata_total.obsm['X_pca'][:, 2], c=colors, s=100)

# Annotate legend with labels
legend_handles = []
for label, color in color_palette.items():
    if label in [ "PC", "immtransB", "CLP", "Ery"]:
        continue
    legend_handles.append(plt.scatter([], [], label=label, color=color))

plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.title('PCA of liCHiC samples')
plt.xlabel(f"PC1 Explained Variability: {adata_total.uns['pca']['variance_ratio'][0]:.3f}")
plt.ylabel(f"PC3 Explained Variability: {adata_total.uns['pca']['variance_ratio'][2]:.3f}")
plt.savefig(f'{outdir}pca_lichic_cell_types_cmp_preb_gcb_pc1_pc3.svg', format='svg', bbox_inches='tight')
plt.clf()

########################################################################################

#Pr clu


#########################################################################################


print("Performing per cluster analysis...") 

adata = ad.read_h5ad("/home/bsc/bsc059153/PROJECTS/BALL/TFM/results/2_clustering_results/clustering_100000_4000_10_raw_all/anndata/100000_4000_10_raw_all_anndata.hdf5")

print(adata.obs['leiden'].unique())

# Ensure the 'leiden' column exists
if 'leiden' in adata.obs.columns:
    # Split the observations into two lists based on the value of 'leiden'
    group1 = adata[adata.obs['leiden'] == '0']
    group2 = adata[adata.obs['leiden'] == '1']
    
else:
    print("The 'leiden' column is not present in the AnnData object.")


group1 = group1.obs_names.tolist()
group2 = group2.obs_names.tolist()

group1_dict = {}
group2_dict = {}

#Create a dictioanry of this

for element in group1:
    group1_dict[element] = 1
for element in group2:
    group2_dict[element] = 1


cells_group1, group1_features = process_demux_file_multi(group1_dict, good_cells, reso, too_close, demux_file)
cells_group2, group2_features = process_demux_file_multi(group2_dict, good_cells, reso, too_close, demux_file)

#####

# Comparison

n_feat_sc, n_feat_lichic, n_feat_inter, similarity_score, common_features_group1 = jaccard_similarity(list(group1_features), list(lichic_relapse.keys()))

correlation, pvalue = spearman_correlation_lichic_neo(cells_relapse, lichic_relapse)

#Print results to stdout

print(f"There are {n_feat_inter} unique contacts from the lichic data ({n_feat_lichic} total unique contacts) that are present in the schic data ({n_feat_sc} different interactions using group1 data)")
print(f"Similarity score: {similarity_score}")
print(f"Spearman correlation: {correlation}")


n_feat_sc, n_feat_lichic, n_feat_inter, similarity_score, common_features_group2 = jaccard_similarity(list(group2_features), list(lichic_relapse.keys()))

correlation, pvalue = spearman_correlation_lichic_neo(cells_diagnose, lichic_relapse)

#Print results to stdout

print(f"There are {n_feat_inter} unique contacts from the lichic data ({n_feat_lichic} total unique contacts) that are present in the schic data ({n_feat_sc} different interactions using group2 data)")
print(f"Similarity score: {similarity_score}")
print(f"Spearman correlation: {correlation}")


# Generate venn diagram

# Get intersection between relapse and diagnose

common_values = set(common_features_group1).intersection(common_features_group2)

# Number of common interactions

print(f"Number of common interactions between sets: {len(list(common_values))}")

# Similarity score BETWEEN sets

similarity_among_sets = len(list(common_values)) / ( len(common_features_group1) + len(common_features_group2) - 2*len(list(common_values)) )

print(f"Similarity between relapse and diagnose: {similarity_among_sets}")

plt.figure(figsize=(10, 6), dpi=100)

# Create a Venn diagram

venn_labels = {'100': 'Common Bins with lichic in group 2', '010': 'Common Bins in with lichic in group 1', '110': 'Intersection'}
venn2(subsets=[set(common_features_group2), set(common_features_group1)], set_labels=('Common Bins in group2', 'Common Bins in group1'))

plt.savefig(f"{outdir}venn_diagram_per_cluster.png")

print("Saving unique contacts into files...")

_, _, group2_feat_in_cell = dict_to_vector(cells_group2)    #Binarize the DipC data and obtain feature dictionary

_, _, group1_feat_in_cell = dict_to_vector(cells_group1)    #Binarize the DipC data and obtain feature dictionary


count_bins_group2 = {k: v for k, v in group2_feat_in_cell.items() if v >= 3}
threshold_bins_group2 = set(count_bins_group2.keys())

count_bins_group1 = {k: v for k, v in group1_feat_in_cell.items() if v >= 3}
threshold_bins_group1 = set(count_bins_group1.keys())


#Here change the group2_features and grouop1_features to include features that are exclusive to each condition. Here we are juist using the whol set of features found in each condition

cells_relapse, relapse_features = process_demux_file_multi(cells_group2, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins = group2_features, threshold_bins = threshold_bins_group2, output_file="../6_gsea/group2_exclusive_features.txt")
cells_relapse, relapse_features = process_demux_file_multi(cells_group1, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins = group1_features, threshold_bins = threshold_bins_group1, output_file="../6_gsea/group1_exclusive_features.txt")

#######################################33

result_dict = {}

for filename in ibed_files:
    # Extracting the prefix before the first underscore
    prefix = filename.split('_')[0]
    if prefix in ["PC", "immtransB", "Ery", "CLP", "relapse1", "diagnose1"]:
        continue
    full_name =  f"../data/lichic_files/{filename}"
    print(full_name)
    result_dict[prefix] = process_lichic_file(full_name, too_close, reso)

lichic_file="/home/bsc/bsc059153/PROJECTS/BALL/TFM/data/lichic_files/BALL_relapse_2_cutoff_5.csv"
lichic_relapse = process_lichic_file(lichic_file, reso, too_close)
result_dict["relapse2"] = lichic_relapse

group2_tags = []
group1_tags = []

for i in adata.obs_names:
    if adata.obs["leiden"][i] == "0":
        group1_tags.append(i)
    else:
        group2_tags.append(i)
        
cells_group1_pseudo, _ =  process_demux_file_multi(group1_dict, good_cells, reso, too_close, demux_file)
        
group1_schic = {}
for key, value in cells_group1_pseudo.items():
    for contact, values in value.items():
        try:
            group1_schic[contact] += values
        except KeyError:
            group1_schic[contact] = values
            
new_group1_schic = {key: value for key, value in relapse_schic.items() if value > 5}

cells_group2_pseudo, _ =  process_demux_file_multi(group2_dict, good_cells, reso, too_close, demux_file)
        
group2_schic = {}
for key, value in cells_group2_pseudo.items():
    for contact, values in value.items():
        try:
            group2_schic[contact] += values
        except KeyError:
            group2_schic[contact] = values
            
new_group2_schic = {key: value for key, value in group2_schic.items() if value > 5}

result_dict["group2"] = new_group2_schic
result_dict["group1"] = new_group1_schic


#########

# Binarize data

""" for key in result_dict:
    for subkey in result_dict[key]:
        # Assign 1 if value > 0, else leave it at 0
        if result_dict[key][subkey] > 1:
            result_dict[key][subkey] = 1
        else:
            result_dict[key][subkey] = 0 """

#########
# Get PCA plot

print("Performing PCA of cell types")

# Extract unique keys from all inner dictionaries
inner_keys = sorted(set().union(*[d.keys() for d in result_dict.values()]))

# Create a list of lists to hold the values
feature_idx = dict((f, i) for i, f in enumerate(inner_keys))
X = sparse.lil_matrix((len(result_dict.keys()), len(inner_keys)), dtype=np.double)
for i, (c, ff) in enumerate(result_dict.items()):
    for f, v in ff.items():
        #print(v)       
        X[i, feature_idx[f]] = int(v)


# Convert dictionary of dictionaries to Anndata object
adata_total = ad.AnnData(X=X.tocsc(), obs=list(result_dict.keys()), 
                  var=np.array(list(inner_keys)), 
                  dtype=X.dtype)

adata_total.var_names  = np.array(list(inner_keys))
adata_total.obs_names = list(result_dict.keys())

del(adata_total.obs[0])
del(adata_total.var[0])

sc.pp.filter_genes(adata_total, min_cells=3, inplace=True)

adata_total.obs['total'] = ""
for cell_type in result_dict.keys():    
    adata_total.obs['total'][cell_type] = sum(result_dict[cell_type].values()) 

adata_total.obs['total'] = adata_total.obs['total'].astype(float)

sc.pp.regress_out(adata_total, keys=['total'])
sc.tl.pca(adata_total, n_comps=10)

######

#Hierearchrial


adata_total.obsm['X_pca'] = adata_total.obsm['X_pca']  # Ensure PCA is in the AnnData object
linkage_matrix = linkage(adata_total.obsm['X_pca'], method='ward')

# Create a dendrogram plot
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=adata_total.obs_names)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Sample Index')
plt.ylabel('Distance')

# Rotate the x-axis labels to be diagonal
plt.xticks(rotation=45, ha='right')
plt.tight_layout()


plt.savefig(f'{outdir}/hierarchical_clustering_dendrogram_with_clusters.png')
plt.show()

######


labels = adata_total.obs.index  # replace 'label' with the actual column name in your data
plt.figure(figsize=(6, 6))

color_palette = {
    'Mon': 'yellow',
    'HSC': 'black',
    'CMP': 'grey',
    'CLP': 'grey',
    'nCD4': 'lightblue',
    'nCD8': 'lightgreen',
    'nB': 'cyan',
    'memB': 'blue',
    'tB': 'mediumblue',
    'PC': 'lime',
    'GCB': 'darkgreen',
    'immtransB': 'aquamarine',
    'PreB': 'orange',
    'PreProB': 'indianred',
    'ProB': 'firebrick',
    'sc_relapse': 'tomato',
    'sc_diagnose': 'pink',
    'group1': 'red',
    'relapse2': 'green',
    'group2': 'purple',
    'Ery': 'olive'
}

# Create a list of colors based on cell types
colors = [color_palette[cell_type] for cell_type in adata_total.obs_names]


# PC1 VS PC2

plt.scatter(adata_total.obsm['X_pca'][:, 0], adata_total.obsm['X_pca'][:, 1], c=colors, s=100)

# Annotate legend with labels
legend_handles = []
for label, color in color_palette.items():
    if label in [ "PC", "immtransB", "CLP", "Ery"]:
        continue
    legend_handles.append(plt.scatter([], [], label=label, color=color))

plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.title('PCA of liCHiC samples')
plt.xlabel(f"PC1 Explained Variability: {adata_total.uns['pca']['variance_ratio'][0]:.3f}")
plt.ylabel(f"PC2 Explained Variability: {adata_total.uns['pca']['variance_ratio'][1]:.3f}")
plt.savefig(f'{outdir}/pca_lichic_cell_types_cluster_pc1_pc2.svg', format='svg', bbox_inches='tight')
plt.clf()

# PC1 VS PC3

plt.scatter(adata_total.obsm['X_pca'][:, 0], adata_total.obsm['X_pca'][:, 2], c=colors, s=100)

# Annotate legend with labels
legend_handles = []
for label, color in color_palette.items():
    if label in [ "PC", "immtransB", "CLP", "Ery"]:
        continue
    legend_handles.append(plt.scatter([], [], label=label, color=color))

plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.title('PCA of liCHiC samples')
plt.xlabel(f"PC1 Explained Variability: {adata_total.uns['pca']['variance_ratio'][0]:.3f}")
plt.ylabel(f"PC3 Explained Variability: {adata_total.uns['pca']['variance_ratio'][2]:.3f}")
plt.savefig(f'{outdir}/pca_lichic_cell_types_cluster_pc1_pc3.svg', format='svg', bbox_inches='tight')
plt.clf()

####################3

""" # Per cell type analysis


cluster_correlations = {}
cluster_correlations["group1"] = {}
cluster_correlations["group2"] = {}


for cell_type in adata_total.obs_names:
    if cell_type in ["group1", "group2"]:
        continue
    else:
        
        cluster_correlations["group1"][cell_type], _ = spearman_correlation_multi_lichic(result_dict["group1"], result_dict[cell_type])
        cluster_correlations["group2"][cell_type], _ = spearman_correlation_multi_lichic(result_dict["group2"], result_dict[cell_type])
        

print(cluster_correlations)
# Prepare data for plotting
group1_correlations = list(cluster_correlations["group1"].values())
group2_correlations = list(cluster_correlations["group2"].values())
cell_types = list(cluster_correlations["group1"].keys())

print(group2_correlations)
print(group1_correlations)
# create the dataframe
df = pd.DataFrame({'group1': group1_correlations, 'group2': group2_correlations}, index=cell_types)


# plot
ax = df.plot(kind='bar', figsize=(6, 4), rot=0, title='Case Comparison', ylabel='Values')
plt.savefig(f'{outdir}/correlation_cell_types_histogram.png') """
