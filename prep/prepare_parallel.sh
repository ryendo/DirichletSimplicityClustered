#!/bin/bash
# Define MAX_JOBS at the beginning of the script
MAX_JOBS=60  # Adjust this value as needed
# Get the absolute path of the current directory
SCRIPT_DIR=$(pwd)
# Define source and target directories
SOURCE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group/Intlab_V12"
TARGET_BASE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group"
# Ensure the target base directory exists
mkdir -p "$TARGET_BASE_DIR"

# ▼▼▼ Change: Path to check for list_j.csv is updated ▼▼▼
# Check if list_j.csv exists
if [ ! -f "${SCRIPT_DIR}/prep/list_j.csv" ]; then
    echo "Error: 'prep/list_j.csv' does not exist. Please run prep.sh first."
    exit 1
fi

# Copy the directory and rename it for each job
for ((i=1; i<=MAX_JOBS; i++)); do
    TARGET_DIR="${TARGET_BASE_DIR}/Intlab_V12_no${i}"
    if [ -d "$TARGET_DIR" ]; then
        rm -rf "$TARGET_DIR"  # Remove existing directory without prompt
    fi
    cp -r "$SOURCE_DIR" "$TARGET_DIR"
done
# (The rest of the script is unchanged)