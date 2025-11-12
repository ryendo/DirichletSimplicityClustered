#!/bin/bash
# This script runs the main computation in a robust, parallel, and restartable manner.

# --- Configuration ---
MAX_JOBS=1  # Adjust this value as needed
FUNC_SCRIPT="func_algo2"

# --- Setup ---
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
mkdir -p "${SCRIPT_DIR}/Each_Process/log"
LOCKFILE="${SCRIPT_DIR}/prep/algo2_list_j.lock"

list_idle_process=()
for ((i=1; i<=MAX_JOBS; i++)); do
    list_idle_process+=($i)
done

# --- Main Execution ---
echo "Launching ${MAX_JOBS} parallel worker processes..."

for n in "${list_idle_process[@]}"; do
    (
        cd "${SCRIPT_DIR}/Each_Process" || exit 1
        export DISPLAY=

        matlab -nodisplay -nosplash -nodesktop -r "try; addpath('${SCRIPT_DIR}/Each_Process'); my_intlab_mode_config(${n}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);" > "log/process_no${n}.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "FATAL: Initial setup (my_intlab_mode_config) failed for worker ${n}. Check log/process_no${n}.log"
            exit 1
        fi

        while true; do
            j=""
            line_number=""
            # --- FIX: Use { ... } group instead of (...) subshell to preserve variable scope ---
            {
                flock -x 200
                j_line=$(awk -F',' 'NF==2 && $2==0 {print NR","$1; exit}' "${SCRIPT_DIR}/prep/algo2_list_j.csv")
                
                if [ -n "$j_line" ]; then
                    line_number=$(echo "$j_line" | cut -d',' -f1)
                    j=$(echo "$j_line" | cut -d',' -f2)   # for example j="1" or j="2"
                    sed -i.bak "${line_number}s/.*/${j},2/" "${SCRIPT_DIR}/prep/algo2_list_j.csv"
                fi
            } 200>"$LOCKFILE"
            
            if [ -z "$j_line" ]; then
                echo "Worker ${n}: No more jobs to process. Exiting."
                break
            fi

            echo "Worker ${n}: Starting job j=${j}..."

            "${SCRIPT_DIR}/prep/run_matlab.sh" "$n" "$FUNC_SCRIPT" "$j" >> "log/process_no${n}.log" 2>&1 &
            wait $!
            matlab_exit_code=$?
        done
        exit 0
    ) &
done

echo "All worker processes launched. Waiting for completion... (You can safely close this terminal)"
wait
echo "All processes have completed."
