#!/bin/bash
#
# This script performs all necessary pre-computation steps for parallel execution.
# 1. Generate the task list (list_j.csv)
# 2. Set up the parallel execution environment (by copying INTLAB)

echo "--- Starting project preparation ---"
echo ""

# Get the directory where this script resides (= project root).
# This ensures that paths work correctly regardless of where the script is called from.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Step 1: Generate the task list (list_j.csv)
echo "[Step 1/2] Generating task list 'list_j.csv'..."
matlab -nodisplay -nosplash -nodesktop -r "try; run('${SCRIPT_DIR}/prep/make_indices_algo2.m'); catch ME; disp(getReport(ME)); exit(1); end; exit(0);"

# ▼▼▼ Change: Path to check for the generated file is updated ▼▼▼
# Check for successful execution
if [ ! -f "${SCRIPT_DIR}/prep/list_j.csv" ]; then
    echo "Error: Failed to generate 'prep/list_j.csv'. Please check the MATLAB output."
    exit 1
fi
echo "=> 'prep/list_j.csv' has been successfully generated in prep/ directory."
echo ""

# Step 2: Set up the parallel execution environment
echo "[Step 2/2] Preparing INTLAB environments for parallel execution..."
bash "${SCRIPT_DIR}/prep/prepare_parallel.sh"

echo ""
echo "--- All preparation steps are complete. ---"