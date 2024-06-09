#!/bin/bash

# Specify the pattern for files
pattern="/gpfs/scratch/bsc08/bsc08786/PROJECTS/HUMAN/BALL/BALL_diagnose/results/H*/"

# Specify the other file
another_file="/gpfs/scratch/bsc08/bsc08786/PROJECTS/HUMAN/BALL/BALL_diagnose/results/IgG_CUTnRUN/alignment/BALL_diagnose_IgG_2_CUTnRUN.filt.nodup.bam"

# Loop over the files matching the pattern
while IFS= read -r -d '' file; do
    # Extract the filename without the directory path
    filename=$(basename "$file")

    # Extract the substrings after the first and second underscores
    diagnose=$(echo "$filename" | cut -d '_' -f2)
    H2AK119ub=$(echo "$filename" | cut -d '_' -f3)

    # Print the desired output
    echo -e "$diagnose\t$H2AK119ub\t$file\t$another_file"
done < <(find /gpfs/scratch/bsc08/bsc08786/PROJECTS/HUMAN/BALL/BALL_diagnose/results/H*/ -name "*_2_*filt.nodup.bam" -print0)
