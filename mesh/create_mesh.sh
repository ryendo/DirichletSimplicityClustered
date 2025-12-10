#!/usr/bin/bash

gmsh_command=$1 
mesh_path=$2
echo "mesh_path", mesh_path
# Force gmsh to use system C++ runtime (Ubuntu) rather than the ones overwritten in matlab env
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

# Optional: clear MATLAB’s library path if called from MATLAB
unset LD_LIBRARY_PATH

${gmsh_command} ${mesh_path}temp.geo -2 -format msh2 -o ${mesh_path}temp.msh 
