import pandas as pd
import sys

# Function to check if two ranges overlap
def ranges_overlap(start1, end1, start2, end2):
    return end1 >= start2 and end2 >= start1

# Ensure the correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python process_data.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the first file containing chromosome ranges
file1 = pd.read_csv("promoter_regions_Homo_sapiens.GRCh38.104_chr.bed", sep="\t", header=0, usecols=[0, 1, 2, 6], names=["chromosome", "range_start", "range_end", "gene"])

# Convert file1 to a dictionary for faster lookup
promoter_dict = {}
for index, row in file1.iterrows():
    chrom = row['chromosome']
    if chrom not in promoter_dict:
        promoter_dict[chrom] = []
    promoter_dict[chrom].append((row['range_start'], row['range_end'], row['gene']))

#print(promoter_dict)

# Read the second file containing Dip-C data
file2 = pd.read_csv(input_file, sep="\t", header=None, names=["info", "chrom1", "pos1", "foo1", "foo2", "read1_start", "read1_end", "chrom2", "pos2", "foo3", "foo4", "read2_start", "read2_end", "ct"])

# List to hold results
results = []

# Iterate through each row in file2
for index2, row2 in file2.iterrows():
    print(index2)
    chrom1 = row2['chrom1']
    chrom2 = row2['chrom2']
    start1 = row2['read1_start'] - 500
    end1 = row2['read1_end'] + 500
    start2 = row2['read2_start'] - 500
    end2 = row2['read2_end'] + 500

    gene_found = False

    if chrom1 in promoter_dict:
        for range_start, range_end, gene in promoter_dict[chrom1]:
            if ranges_overlap(start1, end1, range_start, range_end):
                results.append(row2.tolist() + [gene])
                gene_found = True
                break

    if not gene_found and chrom2 in promoter_dict:
        for range_start, range_end, gene in promoter_dict[chrom2]:
            if ranges_overlap(start2, end2, range_start, range_end):
                results.append(row2.tolist() + [gene])
                break

# Write results to output file
result_df = pd.DataFrame(results, columns=list(file2.columns) + ['gene'])
result_df.to_csv(output_file, sep='\t', index=False)

