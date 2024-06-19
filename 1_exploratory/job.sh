#!/bin/bash
#SBATCH --job-name=translocations
#SBATCH --output=translocations.out
#SBATCH --error=translocations.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00



# Run the Python script with parameters
python translocation_finde.py ../results/2_clustering_results/clustering_1000000_5000_5_raw_all/anndata/1000000_5000_5_raw_all_anndata.hdf5
