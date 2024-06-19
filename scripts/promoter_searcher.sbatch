#!/bin/bash

#SBATCH --job-name="promoter_searcher"
#SBATCH --output=/ijc/LABS/JAVIERRE/DATA/SCHIC/TFM_NICO/promoter_search_dipc/promoter%j.err
#SBATCH --error=/ijc/LABS/JAVIERRE/DATA/SCHIC/TFM_NICO/promoter_search_dipc/promoter%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=1Gb
#SBATCH --qos=hpc

# Get the input and output file names from the command line arguments
input_file=$1
output_file=$2

# Your processing command here, for example:
python promoter_schic_searcher_neo_gene.py "$input_file" "$output_file"

