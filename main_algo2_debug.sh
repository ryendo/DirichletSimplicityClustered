#!/bin/bash
#
# main_algo2.sh (Debug Version)
# This version adds 'echo' statements to trace the script's execution flow.
#

echo "--- Starting main_algo2.sh in DEBUG mode ---"

# --- Fix 1: Use a reliable method to get the script's directory ---
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo "[Debug] Project root directory set to: ${SCRIPT_DIR}"

MAX_JOBS=4
FUNC_SCRIPT="func_algo2"

# --- Fix 2: Add a check before creating the log directory ---
echo "[Debug] Preparing to create log directory: ${SCRIPT_DIR}/Each_Process/log"
mkdir -p "${SCRIPT_DIR}/Each_Process/log"
if [ ! -d "${SCRIPT_DIR}/Each_Process/log" ]; then
    echo "Error: Failed to create log directory. Check permissions."
    exit 1
fi
echo "[Debug] Log directory is ready."

LOCKFILE="${SCRIPT_DIR}/prep/list_j.lock"

# --- Fix 3: Run only ONE worker for the test, not in the background ---
echo "[Debug] Will launch a single test worker (n=1)."
n=1 # Test with only the first worker

# The subshell logic starts here
(
    echo "[Worker ${n}] Sub-process started. PID is $$."

    # Change to the Each_Process directory
    TARGET_CD_DIR="${SCRIPT_DIR}/Each_Process"
    echo "[Worker ${n}] Attempting to change directory to: ${TARGET_CD_DIR}"
    cd "${TARGET_CD_DIR}" || { echo "[Worker ${n}] Error: Failed to change directory. Exiting."; exit 1; }
    echo "[Worker ${n}] Successfully changed directory."

    # Set DISPLAY variable to prevent graphics usage
    export DISPLAY=

    # Execute my_intlab_mode_config(n) in MATLAB
    LOG_FILE="log/process_no${n}.log"
    echo "[Worker ${n}] Running initial setup: my_intlab_mode_config(${n})..."
    echo "[Worker ${n}] Log output will be written to: ${TARGET_CD_DIR}/${LOG_FILE}"

    matlab_output=$(matlab -nodisplay -nosplash -nodesktop -r "try; my_intlab_mode_config(${n}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);" 2>&1)
    matlab_exit_code=$?
    
    echo "[Worker ${n}] Initial setup finished with exit code: ${matlab_exit_code}."

    # Append the output to the log file
    echo "$matlab_output" > "${LOG_FILE}" # Use '>' for the first write to create the file

    if [ $matlab_exit_code -ne 0 ]; then
        echo "[Worker ${n}] Error during initial setup. Check the log file for details."
        echo "Error: my_intlab_mode_config(${n}) failed in process ${n}." >> "${LOG_FILE}"
        echo "MATLAB Error Output:" >> "${LOG_FILE}"
        echo "$matlab_output" >> "${LOG_FILE}"
        exit 1
    fi
    echo "[Worker ${n}] Initial setup successful. Starting main loop..."

    # The rest of the script (while loop) would go here,
    # but for a simple test, we can stop here.
    echo "[Worker ${n}] Debug test finished successfully."

)
# Note: The '&' is removed for this test to see output directly.

echo ""
echo "--- Debug script finished. Check the output above. ---"