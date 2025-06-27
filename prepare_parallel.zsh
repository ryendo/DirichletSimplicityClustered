#!/bin/zsh

# Define MAX_JOBS at the beginning of the script
MAX_JOBS=4  # Adjust this value as needed

# Get the absolute path of the current directory
SCRIPT_DIR=$(pwd)

# Define source and target directories
SOURCE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group/Intlab_V12"
TARGET_BASE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group"

# Ensure the target base directory exists
mkdir -p "$TARGET_BASE_DIR"

# Check if list_j.csv exists
if [ ! -f "${SCRIPT_DIR}/list_j.csv" ]; then
    echo "Error: list_j.csv does not exist."
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

# Initialize list_idle_process with integers from 1 to MAX_JOBS
list_idle_process=()
for ((i=1; i<=MAX_JOBS; i++)); do
    list_idle_process+=($i)
done

# Define the function script name
FUNC_SCRIPT="func_main_bounds"