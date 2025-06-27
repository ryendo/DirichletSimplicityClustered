#!/bin/bash

# MATLAB script to configure intlab
CONFIG_SCRIPT="my_intlab_mode_config"

# Function script to call with j values
FUNC_SCRIPT="func_main_bounds_monotonicity"

# Maximum number of concurrent jobs
MAX_JOBS=80

# Function to run MATLAB with given j values
run_matlab() {
    local j_list=("$@")
    j_list_str=$(IFS=,; echo "${j_list[*]}")
    nohup matlab -nodisplay -nosplash -r "run('$CONFIG_SCRIPT'); $FUNC_SCRIPT([$j_list_str]); exit;" &
}

# Create arrays of j values in chunks of 10 and execute jobs
for ((i=1; i<=19; i+=1)); do
    j_list=()
    for ((j=i; j<i+1 && j<=19; j++)); do
        j_list+=($j)
    done
    
    # Run MATLAB function with the current chunk of j values
    run_matlab "${j_list[@]}"
    
    # Wait if the number of concurrent jobs reaches the limit and check for finished jobs
    while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
        # Remove finished jobs
        jobs -rp | xargs -n 1 -r kill -0 2>/dev/null || wait
    done
done

# Wait for all background jobs to finish
wait
