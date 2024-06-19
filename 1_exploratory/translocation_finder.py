import anndata as ad
import numpy as np
import scanpy as sc
import sys
sys.path.append("../scripts")
from hic_functions import hic_zoom, translocation_contact_finder, hic_zoom_zoom, hic_zoom_fusion
from filter_scHiC import get_cells, cis_ratio, do_plot
from CUTAG_lib.cutag.utilities.clustering import wanted_leiden
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score
import pandas as pd
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import scipy as sp
# Disable warnings
warnings.filterwarnings("ignore")

outdir="../results/1_exploratory"
anndata_object=sys.argv[1]

#################################################################

#Set parameters

reso = 1_000_000
trans_reso = 1_000_000
min_count_per_cell = 5000
too_close = 10_000 
min_counts_feature = 5
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

all_cells = {key: value for key, value in my_dict.items()}

#translocation_file_pd = pd.read_csv('../data/Miltelman_BALL_gene-fusions.csv', delimiter=',', header=0, names=['BALL',  'fusion_protein', 'refs', 'gene1', 'gene2', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])

translocation_file_pd = pd.read_csv('../data/Miltelman_BALL_translocations.csv', delimiter=',', header=0, names=['BALL',  'translocation', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'refs'])

adts = ad.read_h5ad(anndata_object)

print("Getting hic map plot ...")

for index, row in translocation_file_pd.iterrows():
    translocation = row['translocation']
    chr1 = row['chr1']
    chr2 = row['chr2']
    start1 = int(row['start1'])-5000000
    end1 = int(row['end1'])+5000000
    start2 = int(row['start2'])-5000000
    end2 = int(row['end2'])+5000000
    
    print(f"Processing translocation {translocation} between {chr1} and {chr2}...")
    #if translocation!="t(1;5)(p22;p22)":
        #continue
    #cells2 = hic_zoom_fusion(all_cells, good_cells, demux_file, 1_000_000, 100000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, chr1 ,chr2 , outdir)
    cells2, _, chromosome_position= hic_zoom(all_cells, good_cells, demux_file, 1_000_000, 100_000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, outdir)
    cells2, relapse_matrix, chromosome_position = hic_zoom(relapse_cells, good_cells, demux_file, 1_000_000, 100_000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, outdir)
    cells2, diagnose_matrix, chromosome_position = hic_zoom(diagnose_cells, good_cells, demux_file, 1_000_000, 100_000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, outdir)
    
    ###################################
    
    target_sum = np.sum(relapse_matrix)

    # Compute the current sums of the matrices  
    sum1 = np.sum(diagnose_matrix)
    sum2 = np.sum(relapse_matrix)

    # Compute scaling factors
    scale_factor1 = target_sum / sum1
    scale_factor2 = target_sum / sum2

    # Scale the matrices
    normalized_matrix1 = diagnose_matrix * scale_factor1
    normalized_matrix2 = relapse_matrix * scale_factor2
    
    sigma = 2 # gaussian smoothing, the lower the more normal (0 is doing nothing)

    m1 = sp.ndimage.filters.gaussian_filter(normalized_matrix1, sigma)
    m2 = sp.ndimage.filters.gaussian_filter(normalized_matrix2, sigma)
    
    difference_matrix = m2 - m1 
    
    
    # Normalized matrix1
    plt.subplot(3, 1, 1)  # Change to (1, 3, 1)
    plt.imshow(normalized_matrix1, vmax=10)
    plt.title('Normalized Diagnose')
    plt.colorbar()
    plt.xticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)
    plt.yticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)

    # Normalized matrix2
    plt.subplot(3, 1, 2)  # Change to (1, 3, 2)
    plt.imshow(normalized_matrix2, vmax=10)
    plt.title('Normalized Relapse') 
    plt.colorbar()
    plt.xticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)
    plt.yticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)

    # Difference matrix
    plt.subplot(3, 1, 3)  # Change to (1, 3, 3)
    plt.imshow(difference_matrix, cmap='coolwarm', vmax=10, vmin=-10)
    plt.title('Difference Matrix')
    plt.colorbar()
    plt.xticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)
    plt.yticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)

    plt.tight_layout()

    # Save the figure
    plt.savefig(f'{outdir}/{translocation}_all_matrix_plots.png', dpi=300)
    plt.clf()
    ###################################
    
    trans_reso=1_000_000
    translocation_contact_finder(adts, chr1, start1, end1, chr2, start2, end2, trans_reso, cells2) 

for cell, features in all_cells.items():
    try:
        adts.obs.loc[cell, "BALL_translocations"] = adts.obs.loc[cell, "BALL_translocations"] / adts.obs.loc[cell, "trans"]
    except KeyError:
        continue
    
sc.pl.umap(adts, color=['group','BALL_translocations'], size=200, save=f"_BALL_translocation_embedding.svg")


true_labels = adts.obs["leiden"]  # Replace with your actual class labels

adts.obs['BALL_translocations_binary'] = ((adts.obs['BALL_translocations'] > adts.obs['BALL_translocations'].median())).astype(int)
variable_values = adts.obs['BALL_translocations_binary'] 

nmi = normalized_mutual_info_score(true_labels, variable_values)
print(f"Normalized Mutual Information (leiden) for 'BALL_translocations':", nmi)

true_labels = adts.obs['group_binary'] = np.where(adts.obs["group"] == "diagnose", 1, 0)
nmi = normalized_mutual_info_score(true_labels, variable_values)
print(f"Normalized Mutual Information (group) for 'BALL_translocations':", nmi)

#################################################################

data = adts.obs[['leiden', 'BALL_translocations']]

# Convert to a pandas DataFrame (optional, but can be convenient)
df = pd.DataFrame(data)

print(df)

# Plot the boxplot
plt.figure(figsize=(6, 10))
sns.boxplot(x='leiden', y='BALL_translocations', data=df)
plt.title('BALL Translocations by Cluster')
plt.savefig(f'{outdir}/ball_translocations_by_cluster.svg')
plt.show()

# Statistical test (t-test in this case, assuming normal distribution)
group1 = df[df['leiden'] == '0']['BALL_translocations']
group2 = df[df['leiden'] == '1']['BALL_translocations']

print(group1)
print(group2)

t_stat, p_value = ttest_ind(group1, group2)

print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

# Check significance
if p_value < 0.05:
    print("The difference between groups is statistically significant.")
else:
    print("The difference between groups is not statistically significant.")

#################################################################


""" adts.obs["BALL_translocations"] = 0
translocation_file_pd = pd.read_csv('../data/Miltelman_BALL_gene-fusions.csv', delimiter=',', header=0, names=['BALL',  'fusion_protein', 'refs', 'gene1', 'gene2', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
for index, row in translocation_file_pd.iterrows():
    translocation = row['fusion_protein']
    chr1 = row['chr1']
    chr2 = row['chr2']
    print(f"Processing translocation {translocation} between {chr1} and {chr2}...")
    #if translocation!="t(1;5)(p22;p22)":
        #continue
    cells2 = hic_zoom_fusion(all_cells, good_cells, demux_file, 1_000_000, 100000, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, chr1 ,chr2 , outdir)
    trans_reso=1_000_000
    translocation_contact_finder(adts, chr1, start1, end1, chr2, start2, end2, trans_reso, cells2) 

print("Getting UMAP with translocation...")

for cell, features in all_cells.items():
    try:
        adts.obs.loc[cell, "BALL_translocations"] = adts.obs.loc[cell, "BALL_translocations"] / adts.obs.loc[cell, "trans"]
    except KeyError:
        continue
        
        
sc.pl.umap(adts, color=['group','BALL_translocations'], size=200, save=f"_BALL_gene_fusion_embedding.svg")

true_labels = adts.obs["leiden"]  # Replace with your actual class labels

adts.obs['BALL_translocations_binary'] = ((adts.obs['BALL_translocations'] > adts.obs['BALL_translocations'].median())).astype(int)
variable_values = adts.obs['BALL_translocations_binary'] 

nmi = normalized_mutual_info_score(true_labels, variable_values)
print(f"Normalized Mutual Information (leiden) for 'BALL_translocations':", nmi)

true_labels = adts.obs['group_binary'] = np.where(adts.obs["group"] == "diagnose", 1, 0)
nmi = normalized_mutual_info_score(true_labels, variable_values)
print(f"Normalized Mutual Information (group) for 'BALL_translocations':", nmi) """
 
###############################################################

# BOXPLOTSS PER CONDITION AND PER TRNALSOCATION + THE SUM OF ALL