#!/bin/bash

# Script to sync a file to the remote server via rsync with SSH jump.
# It transfers files into the shoe/build/Reconstruction dir in tier1 (where the macros are).
# Usage: ./sync_file.sh /path/to/local/file

# Check if a file path is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/local/file"
  exit 1
fi

# Get the source file path from the argument
SOURCE="$1"

# Define the destination and bastion
DEST="lpierfederici@ui01-foot.cr.cnaf.infn.it:/opt/exp_software/foot/lpierfederici/shoe/shoe_24.0/shoe/build/Reconstruction"
BASTION="lpierfederici@bastion.cnaf.infn.it"

# Execute rsync with SSH jump
rsync -azvh --progress -e "ssh -J $BASTION" "$SOURCE" "$DEST"
