'''

This file contains some general functions used in Dip-C analysis

'''

from collections import Counter
from scipy.stats import spearmanr
from scipy import sparse
import numpy as np
import anndata as ad
from matplotlib import pyplot as plt

# FUNCTIONS #################################################################################################

def process_demux_file_multi(cells_subset, good_cells, reso, too_close, demux_file, use_trans=True, lichic_bins=None, threshold_bins = None, output_file=None):
    
    if use_trans:
        check_trans = lambda x1, x2: False
    else:
        check_trans = lambda x1, x2: x1 != x2
        
    cells_new = {}
    fh = open(demux_file)
    all_features_total = set()
    trans_reso = reso * 10
    
    cchar = 0
    for line in fh:
        if line.startswith('#'):
            cchar += len(line)
            continue
        break
    fh.seek(cchar)
    
    output_fh = None
    if output_file:
        output_fh = open(output_file, 'w')
        
    for line in fh:
        _, c1, b1, s1, l1, rs1, re1, c2, b2, s2, l2, rs2, re2, ct, gene = line.split()
        # id, crom1, pos1, strand1, length1, start1, end1, crom2, pos2, strand2, length2, start2, end2, cell-tag

        # Filter reads
        if ct not in good_cells.keys():
            continue

        if ct not in cells_subset.keys():
            continue

        if check_trans(c1, c2):
            continue

        if c1 == "chrY" or c2 == "chrY":
            continue
        
        b1, l1, b2, l2, rs1, re1, rs2, re2 = map(float, [b1, l1, b2, l2, rs1, re1, rs2, re2])

        # Define exact position
        p1 = b1 + l1 / 2
        p2 = b2 + l2 / 2

        # Remove diagonal and define features
        if c1 == c2:
            if abs(int(p1) - int(p2)) <= too_close:
                continue

            p1r1 = int(int(rs1) // reso)
            p2r1 = int(int(re1) // reso)
            p1r2 = int(int(rs2) // reso)
            p2r2 = int(int(re2) // reso)

            if p1r1 == p1r2 or p1r1 == p2r2 or p2r1 == p1r2 or p2r1 == p2r2:
                continue
        else:
            p1r1 = int(int(rs1) // trans_reso)
            p2r1 = int(int(re1) // trans_reso)
            p1r2 = int(int(rs2) // trans_reso)
            p2r2 = int(int(re2) // trans_reso)

        feature1 = f"{c1}:{p1r1}_{c2}:{p1r2}"
        feature2 = f"{c1}:{p2r1}_{c2}:{p2r2}"
        feature3 = f"{c1}:{p1r1}_{c2}:{p2r2}"
        feature4 = f"{c1}:{p2r1}_{c2}:{p1r2}"

        iteration_features = set()
        iteration_features.add(feature1)
        iteration_features.add(feature2)
        iteration_features.add(feature3)
        iteration_features.add(feature4)

        all_features_total.add(feature1)
        all_features_total.add(feature2)
        all_features_total.add(feature3)
        all_features_total.add(feature4)

        # Fill dictionary
        for feature in iteration_features: 
            try:
                cells_new[ct][feature] += 1
            except KeyError:
                try:
                    cells_new[ct][feature] = 1
                except KeyError:
                    cells_new[ct] = {feature: 1}
    
        # Here only print counts that are present in lichic bins and appear in at leats x amount of cells specified in the previous line (3)
        if lichic_bins and (feature1 in lichic_bins or feature2 in lichic_bins or feature3 in lichic_bins or feature4 in lichic_bins):
            
            if threshold_bins and (feature1 in threshold_bins or feature2 in threshold_bins or feature3 in threshold_bins or feature4 in threshold_bins):
                
                if output_fh:
                    #print("Printing contact to file...")
                    output_fh.write(f"{gene}\n")

    fh.close()
    if output_fh:
        output_fh.close()
    
    return cells_new, all_features_total




################################################################################################################

def process_lichic_file(lichic_file, too_close, reso):
    
    fh = open(lichic_file)
    lichic_contact_bins = []
    trans_reso = reso*10
    
    for line in fh:
        if line.startswith('bait_chr'):
            continue
        c1, bait_start, bait_end, _, c2, other_start, other_end, _, _, _ = line.strip().split('\t')
    
        # Get the middle of the read
    
        pc1 = (int(bait_start) + int(bait_end)) // 2
        pc2 = (int(other_start) + int(other_end)) // 2
    
    
        if c1 == c2:
        
            if abs(int(pc1) - int(pc2)) <= too_close:
                continue

            p1 = int(int(pc1) // reso)
            p2 = int(int(pc2) // reso)
        
            if p1 == p2:
                continue
        
        else:
        
            p1 = int(int(pc1) // trans_reso)
            p2 = int(int(pc2) // trans_reso)
        
    
        if c1 > c2 or (c1 == "X" and c2 != "X"):
        
            feature = f"chr{c2}:{p2}_chr{c1}:{p1}"
            
        elif c1 == c2:
        
            if p1 < p2:
                feature = f"chr{c1}:{p1}_chr{c2}:{p2}"
            else:
                feature = f"chr{c1}:{p2}_chr{c2}:{p1}"
            
        else:
        
            feature = f"chr{c1}:{p1}_chr{c2}:{p2}"
        

        lichic_contact_bins.append(feature)
        
    lichic_dict = Counter(lichic_contact_bins)
    
    return lichic_dict

################################################################################################################


def get_adts_object(n_cells, feat):
    
    '''
    Given a cell dictionary and a list of features it generates the sparse matrix and the anndata object
    '''

    feature_idx = dict((f, i) for i, f in enumerate(feat))
    X = sparse.lil_matrix((len(n_cells), len(feat)), dtype=np.double)
    for i, (c, ff) in enumerate(n_cells.items()):

        for f, v in ff.items():
            X[i, feature_idx[f]] = v
            
    adata = ad.AnnData(X=X.tocsc(), obs=list(n_cells.keys()), 
        var=np.array(list(feature_idx.keys())), 
        dtype=X.dtype)
    #print(adts)
    adata.var_names = np.array(list(feature_idx.keys()))
    adata.obs_names = list(n_cells.keys())
    del(adata.obs[0])
    del(adata.var[0])
    
    return adata

#################################################################################################################

def jaccard_similarity(list1, list2):
    
    '''
    Given two list of contact features it returns the Jaccard similarity between the two sets
    '''
    set1 = set(list1)
    set2 = set(list2)
    intersection = set1.intersection(set2)
    union = len(set1.union(set2))
    similarity = len(intersection) / union
    return len(set1), len(set2), len(intersection), similarity, intersection

################################################################################################################

def dict_to_vector(dictionary):
    
    dict_contacts = {} 
    
    for dict_cells, dict_contv in dictionary.items():
        for key, value in dict_contv.items() :
            #If we have already set the key before, but find a value thats higher, assign that new value (this method we always take the highest value found for a feature)
            #if key in dict_contacts:
            #   if value > dict_contacts[key]:
                    #dict_contacts[key] = value
                #else:
                    #continue
            if key in dict_contacts:   
                dict_contacts[key] += 1
            else:
                dict_contacts[key] = 1
            #else:
                #dict_contacts[key] = value
                
    keys = list(dict_contacts.keys())
    values = list(dict_contacts.values())
    return keys, values, dict_contacts

################################################################################################################


def spearman_correlation_lichic_neo(dict1, dict2):     #Compare cell clusters to lichic dictionary
    
    _, _, dict_contact = dict_to_vector(dict1)


    
    # Get the common keys present in both dictionaries
    common_keys = set(dict_contact.keys()) & set(dict2.keys())
    print(len(common_keys))
    # Filter dictionaries to include only common keys
    dict1_filtered = {k: dict_contact[k] for k in common_keys}
    dict2_filtered = {k: dict2[k] for k in common_keys}
    
    # Extract values from the filtered dictionaries
    vector1 = list(dict1_filtered.values())
    vector2 = list(dict2_filtered.values())
    
    # Calculate Spearman correlation
    correlation, p_value = spearmanr(vector1, vector2)
    
    return correlation, p_value


#################################################################################################################

def spearman_correlation_multi_lichic_neo(dict1, dict2):  #Compare lichic dictionaries using common keys
    
    # Get the common keys present in both dictionaries
    common_keys = set(dict1.keys()) & set(dict2.keys())
    #print(len(common_keys))
    # Filter dictionaries to include only common keys
    dict1_filtered = {k: dict1[k] for k in common_keys}
    dict2_filtered = {k: dict2[k] for k in common_keys}
    
    # Extract values from the filtered dictionaries
    vector1 = list(dict1_filtered.values())
    vector2 = list(dict2_filtered.values())
    print(vector1)
    # Calculate Spearman correlation
    correlation, p_value = spearmanr(vector1, vector2)
    
    return correlation, p_value

#################################################################################################################

def spearman_correlation_multi_lichic(dict1, dict2):    #Compare lichic dictionaries between one another
    #print("Comparing...")
    keys1 = list(dict1.keys())
    values1 = list(dict1.values())


    key_set = set(keys1)
    vector1 = [values1[keys1.index(k)] if k in key_set else 0 for k in dict2.keys()]   # Here: minimize 0s, if the contact is not present in lichic dont use it
    vector2 = list(dict2.values())


    correlation, p_value = spearmanr(vector1, vector2)
    return correlation, p_value

#################################################################################################################

def jaccard_similarity_boot(list1, list2):
    intersection = 0
    set2 = set(list2)
    
    for element in list1:
        if element in set2:
            intersection += 1
            
    union = len(list1) + len(set2) - intersection
    similarity = intersection / union
    
    return similarity


##################################################################################################################

def flatten_dict(d, parent_key=''):
    items = []
    for k, v in d.items():
        new_key = (parent_key, k) if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key).items())
        else:
            items.append((new_key, v))
    return dict(items)

##################################################################################################################

def hic_zoom(cell_dict, good_cells, demux_file, reso, plot_reso, too_close, biases_cis, biases_trans, translocation_file_pd, translocation, outdir):
    
    if translocation not in set(translocation_file_pd['translocation']):
        print("Please provide a valid translocation")
        return
    
    trans_reso = reso
    chr1 = f"chr{translocation.split(')')[0].split('(')[1].split(';')[0]}"
    chr2 = f"chr{translocation.split(')')[0].split('(')[1].split(';')[1]}"
    chr_list=[chr1, chr2]
       
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
        if line.startswith('bait_chr'):
            continue

        _, c1, b1, s1, l1, _, _, c2, b2, s2, l2, _, _, ct = line.split()
        #c1, bait_start, bait_end, _, c2, other_start, other_end, _, _, _ = line.strip().split('\t')

        if ct not in good_cells:
            continue
        if ct not in cell_dict.keys():
            continue
        
        b1, l1, b2, l2 = map(int, [b1, l1, b2, l2])
        # define exact position
        p1 = b1 + l1 / 2
        #print(p1)
        #print(p2)
        p2 = b2 + l2 / 2
    
        if c1 not in chr_list or c2 not in chr_list:
            continue
    
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

    ######################################################
    
    cells2 = cells
               
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
    
    ######################################################

    feature_idx = dict((f, i) for i, f in enumerate(all_features))
    X = sparse.lil_matrix((len(cells), len(all_features)), dtype=np.double)
    for i, (c, ff) in enumerate(cells.items()):

        for f, v in ff.items():
            X[i, feature_idx[f]] = v
            
    ######################################################


    adts_next = ad.AnnData(X=X.tocsc(), obs=list(cells.keys()), 
                  var=np.array(list(feature_idx.keys())), 
                  dtype=X.dtype)

    adts_next.var_names = np.array(list(feature_idx.keys()))
    adts_next.obs_names = list(cells.keys())
    del(adts_next.obs[0])
    del(adts_next.var[0])
    #print(adts_next)
    
    ########################################################
    
    
    fh = open(demux_file)
    chrom_sizes = {}
    for line in fh:
        if not line.startswith('#'):
            break
        *_, c, v = line.split()
        if c not in chr_list: 
            continue
        chrom_sizes[c] = int(v)
        
    ##########################################################

    chromosomes = dict((k, chrom_sizes[k]//reso + 1) for k in chr_list)

    genome_size = sum(chromosomes.values())

    coord_to_bin = {}

    chromosome_position = {}

    total = 0

    for c,v in chromosomes.items():
        for p in range(0,v):
            coord_to_bin[c,p] = p + total
        chromosome_position[c] = total
        total += v 
    chromosome_position

    chrom_sizes_ratio = {}
    for key, value in chrom_sizes.items():

        if key == "chrMT" or key=="chrY":
            continue
        chrom_sizes_ratio[key] = value//reso +1
    chrom_sizes_ratio



    #######################################################

    vals = adts_next.X.sum(axis=0).tolist()[0]
    
    ratio = trans_reso // reso  
    matrix = np.zeros((genome_size,genome_size))

    for n,a  in enumerate(adts_next.var_names):
        p1, p2 = a.split("_")
        #print(p1,p2)
        p1_c, p1_b = p1.split(":")
        p2_c, p2_b = p2.split(":")

        if p1_c == p2_c:
        
            p1 = coord_to_bin[p1_c, int(p1_b)]
            p2 = coord_to_bin[p2_c, int(p2_b)]

            matrix[p1, p2] = vals[n]
            matrix[p2, p1] = vals[n]

        else:


            new_b1 = (int(p1_b))*ratio
            new_b2 = (int(p2_b))*ratio

            for i in range(ratio):

                for j in range(ratio):

                    try: 
                    
                        matrix[coord_to_bin[p1_c, new_b1 + i], coord_to_bin[p2_c, new_b2 + j]] = vals[n]
                        matrix[coord_to_bin[p2_c, new_b2 + i], coord_to_bin[p1_c, new_b1 + j]] = vals[n]
                    
                    except KeyError: 
                    
                        continue
                    
    plt.figure(figsize=(14, 13))
    plt.imshow(matrix, vmax=5)
    plt.colorbar()
    
    translocation_file_pd.loc[translocation_file_pd['translocation'] == translocation]
     

    x_coord1 = coord_to_bin[chr1, int(translocation_file_pd.loc[translocation_file_pd['translocation'] == translocation]["start1"]) // reso]
    plt.axvline(x_coord1, color="red")
    x_coord2 = coord_to_bin[chr1, int(translocation_file_pd.loc[translocation_file_pd['translocation'] == translocation]["end1"]) // reso]
    plt.axvline(x_coord2, color="red")
        

        
    y_coord1 = coord_to_bin[chr2, int(translocation_file_pd.loc[translocation_file_pd['translocation'] == translocation]["start2"]) // reso]
    plt.axhline(y_coord1, color="red")
    y_coord2 = coord_to_bin[chr2, int(translocation_file_pd.loc[translocation_file_pd['translocation'] == translocation]["end2"]) // reso]
    plt.axhline(y_coord2, color="red")
        

 
    # Set the tick positions and labels
    plt.xticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)
    plt.yticks([i for i in chromosome_position.values()], [i for i in chromosome_position.keys()], rotation=45)

    plt.grid()
    
    plt.savefig(f"{outdir}/{translocation}_hic_map.png")
    
    return cells2


##################################################################################################################


def translocation_contact_finder(adts_object, tchr1, tpos1_start, tpos1_end, tchr2, tpos2_start, tpos2_end, translocation, trans_reso, cells):
    """
    Populate a new observation column with contact counts for the specified contact pair.

    Contact must be a string like chr4_chr17, numerically ordered (X < Y < MT)
    
    Values are in theory raw values
    
    translocations bins must be something like: 1-
    """

    # Initialize a new observation column for the specified contact pair
    adts_object.obs[translocation] = 0
    
    tpos1_start = tpos1_start // trans_reso
    tpos1_end = tpos1_end // trans_reso
    
    tpos2_start = tpos2_start // trans_reso
    tpos2_end = tpos2_end // trans_reso  
    
    
    # Iterate over cells and update contact counts
    for cell, features in cells.items():
        if cell not in adts_object.obs_names:
            continue
        for feature, value in features.items():
            if value < 2:
                continue
            else:
                cbin1, cbin2 = feature.split('_')[0], feature.split('_')[1]
                
                chr1, bin1 = cbin1.split(':')
                
                chr2, bin2 = cbin2.split(':')

                if chr1 != tchr1 or chr2 != tchr2:
                    continue
                else:
                    if int(bin1) < tpos1_end and int(bin1) > tpos1_start and int(bin2) < tpos2_end and int(bin2) > tpos2_start:
                        adts_object.obs.loc[cell, translocation] += value
                        print(feature, value)

                        
                    
        adts_object.obs.loc[cell, translocation] = adts_object.obs.loc[cell, translocation] / adts_object.obs.loc[cell, "trans"]
    #print(adts_object.obs[translocation])
    return adts_object