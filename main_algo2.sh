#!/bin/bash
# Define MAX_JOBS at the beginning of the script
MAX_JOBS=4  # Adjust this value as needed
# Get the absolute path of the current directory
SCRIPT_DIR=$(pwd)
# (Directory definitions are unchanged)
# ...
# Define the function script name
FUNC_SCRIPT="func_algo2"
# Ensure the log directory exists
mkdir -p "${SCRIPT_DIR}/Each_Process/log"

# ▼▼▼ Change 1: Update lock file path to prep/ ▼▼▼
# Define the lock file
LOCKFILE="${SCRIPT_DIR}/prep/list_j.lock"

# Launch processes in parallel
for n in "${list_idle_process[@]}"; do
    (
        # (Setup part is unchanged)
        # ...
        # Start processing j values from list_j.csv
        while true; do
            # Acquire lock
            (
                flock -x 200
                # ▼▼▼ Change 2: Update path for awk to read list_j.csv ▼▼▼
                j_line=$(awk -F',' '$2 != 1 && $2 != 2 {print NR","$1; exit}' "${SCRIPT_DIR}/prep/list_j.csv")
                
                if [ -z "$j_line" ]; then
                    exit 1 # No more j to process
                fi
                
                line_number=$(echo "$j_line" | cut -d',' -f1)
                j=$(echo "$j_line" | cut -d',' -f2)

                # ▼▼▼ Change 3: Update paths for sed to modify list_j.csv ▼▼▼
                sed "${line_number}s/$/,2/" "${SCRIPT_DIR}/prep/list_j.csv" > "${SCRIPT_DIR}/prep/list_j.csv.tmp"
                mv "${SCRIPT_DIR}/prep/list_j.csv.tmp" "${SCRIPT_DIR}/prep/list_j.csv"
            ) 200>"$LOCKFILE"
            
            # Execute MATLAB script (path already corrected in previous step)
            "${SCRIPT_DIR}/prep/run_matlab.sh" $n &
            matlab_exit_code=$?
            
            echo "$matlab_output" >> "log/process_no${n}.log"
            
            # Acquire lock again to update the status
            (
                flock -x 200
                if [ $matlab_exit_code -eq 0 ]; then
                    # ▼▼▼ Change 4: Update paths for sed on success ▼▼▼
                    sed "${line_number}s/,2/,1/" "${SCRIPT_DIR}/prep/list_j.csv" > "${SCRIPT_DIR}/prep/list_j.csv.tmp"
                    mv "${SCRIPT_DIR}/prep/list_j.csv.tmp" "${SCRIPT_DIR}/prep/list_j.csv"
                else
                    echo "Error: ${FUNC_SCRIPT}(${j}) failed in process ${n}." >> "log/process_no${n}.log"
                    echo "MATLAB Error Output:" >> "log/process_no${n}.log"
                    echo "$matlab_output" >> "log/process_no${n}.log"
                    # ▼▼▼ Change 5: Update paths for sed on error ▼▼▼
                    sed "${line_number}s/,2$//" "${SCRIPT_DIR}/prep/list_j.csv" > "${SCRIPT_DIR}/prep/list_j.csv.tmp"
                    mv "${SCRIPT_DIR}/prep/list_j.csv.tmp" "${SCRIPT_DIR}/prep/list_j.csv"
                    exit 1
                fi
            ) 200>"$LOCKFILE"
        done
        exit 0
    ) &
done
# Wait for all background processes to finish
wait