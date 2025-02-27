#!/bin/bash

# Define a dictionary mapping energy levels to their corresponding runs
declare -A energy_runs
energy_runs["100MeV"]="AnaFOOT_TW_Decoded_HIT2022_4767.root AnaFOOT_TW_Decoded_HIT2022_4769.root"
energy_runs["140MeV"]="AnaFOOT_TW_Decoded_HIT2022_4805.root AnaFOOT_TW_Decoded_HIT2022_4806.root"
energy_runs["200MeV"]="AnaFOOT_TW_Decoded_HIT2022_4747.root AnaFOOT_TW_Decoded_HIT2022_4748.root AnaFOOT_TW_Decoded_HIT2022_4749.root"
energy_runs["220MeV"]="AnaFOOT_TW_Decoded_HIT2022_4838.root AnaFOOT_TW_Decoded_HIT2022_4840.root"

# Define the output file prefix
output_prefix="AnaFOOT_TW_Decoded_HIT2022_fragm"

# Loop through the dictionary
for energy in "${!energy_runs[@]}"; do
    # Define the output file name
    output_file="${output_prefix}_${energy}.root"
    
    # Get the list of input files for this energy
    input_files=${energy_runs[$energy]}
    
    # Run hadd -f to merge files
    echo "Merging files for ${energy} into ${output_file}..."
    hadd -f $output_file $input_files
    
    # Check if hadd was successful
    if [ $? -eq 0 ]; then
        echo "Merged successfully into ${output_file}."
    else
        echo "Error merging files for ${energy}."
    fi
done