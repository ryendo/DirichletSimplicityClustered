#!/bin/bash
# This script executes a single MATLAB task.
# It takes 3 arguments:
# $1: The worker number (n) for my_intlab_mode_config
# $2: The MATLAB function name to run (e.g., "func_algo2")
# $3: The argument for the MATLAB function (e.g., a "j" value)

# Check if all 3 arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <worker_n> <func_script> <j_value>"
    exit 1
fi

n=$1
func_script=$2
j=$3

# Construct the MATLAB command using the provided arguments
matlab_command="try; my_intlab_mode_config(${n}); ${func_script}(${j}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);"

# Execute the command, do not just echo it
matlab -nodisplay -nosplash -nodesktop -r "${matlab_command}"