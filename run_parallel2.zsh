#!/bin/zsh

# Define MAX_JOBS at the beginning of the script
MAX_JOBS=4  # Adjust this value as needed

# Get the absolute path of the current directory
SCRIPT_DIR=$(pwd)

# Define source and target directories
SOURCE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group/Intlab_V12"
TARGET_BASE_DIR="${SCRIPT_DIR}/Each_Process/Intlab_Group"

# Initialize list_idle_process with integers from 1 to MAX_JOBS
list_idle_process=()
for ((i=1; i<=MAX_JOBS; i++)); do
    list_idle_process+=($i)
done

# Define the function script name
FUNC_SCRIPT="func_main_bounds"

# Ensure the log directory exists
mkdir -p "${SCRIPT_DIR}/Each_Process/log"

# Define the lock file
LOCKFILE="${SCRIPT_DIR}/list_j.lock"

# Launch processes in parallel
for n in "${list_idle_process[@]}"; do
    (
        # Change to the Each_Process directory
        cd "${SCRIPT_DIR}/Each_Process" || exit 1

        # Set DISPLAY variable to prevent graphics usage
        export DISPLAY=

        # Execute my_intlab_mode_config(n) in MATLAB without graphics
        # Capture the output and error messages
        matlab_output=$(matlab -nodisplay -nosplash -nodesktop -r "try; my_intlab_mode_config(${n}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);" 2>&1)
        matlab_exit_code=$?

        # Append the output to the log file
        echo "$matlab_output" >> "log/process_no${n}.log"

        if [ $matlab_exit_code -ne 0 ]; then
            echo "Error: my_intlab_mode_config(${n}) failed in process ${n}." >> "log/process_no${n}.log"
            echo "MATLAB Error Output:" >> "log/process_no${n}.log"
            echo "$matlab_output" >> "log/process_no${n}.log"
            exit 1
        fi

        # Start processing j values from list_j.csv
        while true; do
            # Acquire lock
            (
                flock -x 200

                # Read list_j.csv line by line to find an unprocessed j
                j_line=$(awk -F',' '$2 != 1 && $2 != 2 {print NR","$1; exit}' "${SCRIPT_DIR}/list_j.csv")
                if [ -z "$j_line" ]; then
                    # No more j to process; exit loop
                    exit 1
                fi

                # Extract line number and j value
                line_number=$(echo "$j_line" | cut -d',' -f1)
                j=$(echo "$j_line" | cut -d',' -f2)

                # Mark j as being processed (status 2)
                sed "${line_number}s/$/,2/" "${SCRIPT_DIR}/list_j.csv" > "${SCRIPT_DIR}/list_j.csv.tmp"
                mv "${SCRIPT_DIR}/list_j.csv.tmp" "${SCRIPT_DIR}/list_j.csv"

            ) 200>"$LOCKFILE"

            # Execute func_main_bounds(j) in MATLAB without graphics
            # Capture the output and error messages
            # matlab_output=$(matlab -nodisplay -nosplash -nodesktop -r "try; my_intlab_mode_config(${n}); ${FUNC_SCRIPT}(${j}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);" 2>&1)
            ./run_matlab.sh $n &
            matlab_exit_code=$?

            # Append the output to the log file
            echo "$matlab_output" >> "log/process_no${n}.log"

            # Acquire lock again to update the status
            (
                flock -x 200

                if [ $matlab_exit_code -eq 0 ]; then
                    # On success, mark j as processed (status 1)
                    sed "${line_number}s/,2/,1/" "${SCRIPT_DIR}/list_j.csv" > "${SCRIPT_DIR}/list_j.csv.tmp"
                    mv "${SCRIPT_DIR}/list_j.csv.tmp" "${SCRIPT_DIR}/list_j.csv"
                else
                    # On error, log the error and reset the status
                    echo "Error: ${FUNC_SCRIPT}(${j}) failed in process ${n}." >> "log/process_no${n}.log"
                    echo "MATLAB Error Output:" >> "log/process_no${n}.log"
                    echo "$matlab_output" >> "log/process_no${n}.log"
                    # Reset the status to unprocessed
                    sed "${line_number}s/,2$//" "${SCRIPT_DIR}/list_j.csv" > "${SCRIPT_DIR}/list_j.csv.tmp"
                    mv "${SCRIPT_DIR}/list_j.csv.tmp" "${SCRIPT_DIR}/list_j.csv"
                    # Exit the process
                    exit 1
                fi

            ) 200>"$LOCKFILE"

        done

        # Exit the process when all j values are processed
        exit 0
    ) &
done

# Wait for all background processes to finish
wait
