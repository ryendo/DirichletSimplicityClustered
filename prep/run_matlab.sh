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

# 1. Get the absolute path of the main project directory (one level up from 'prep')
PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &>/dev/null && pwd)

# 2. Add the path to 'Each_Process' to the MATLAB command
addpath_command="addpath('${PROJECT_DIR}/Each_Process');"
# --- FIX ENDS HERE ---


# Construct the MATLAB command using the provided arguments, now with the addpath command
matlab_command="try; ${addpath_command} my_intlab_mode_config(${n}); ${func_script}(${j}); catch ME; disp(getReport(ME, 'extended')); exit(1); end; exit(0);"

# Execute the command
matlab -nodisplay -nosplash -nodesktop -r "${matlab_command}"
