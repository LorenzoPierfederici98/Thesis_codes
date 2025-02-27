#!/bin/bash

# Define a dictionary mapping energy levels to their corresponding runs
declare -A energy_runs
energy_runs["100MeV"]="AnaFOOT_TW_Decoded_HIT2022_4766.root"
energy_runs["140MeV"]="AnaFOOT_TW_Decoded_HIT2022_4727.root AnaFOOT_TW_Decoded_HIT2022_4801.root AnaFOOT_TW_Decoded_HIT2022_4802.root"
energy_runs["180MeV"]="AnaFOOT_TW_Decoded_HIT2022_4725.root AnaFOOT_TW_Decoded_HIT2022_4726.root"
energy_runs["200MeV"]="AnaFOOT_TW_Decoded_HIT2022_4742.root AnaFOOT_TW_Decoded_HIT2022_4743.root AnaFOOT_TW_Decoded_HIT2022_4744.root AnaFOOT_TW_Decoded_HIT2022_4745.root"
energy_runs["220MeV"]="AnaFOOT_TW_Decoded_HIT2022_4828.root AnaFOOT_TW_Decoded_HIT2022_4830.root AnaFOOT_TW_Decoded_HIT2022_4832.root AnaFOOT_TW_Decoded_HIT2022_4833.root"

# Define the output file prefix
output_prefix="AnaFOOT_TW_Decoded_HIT2022"

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