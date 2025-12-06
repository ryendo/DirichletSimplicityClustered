#!/usr/bin/bash

mesh_path=$1

# Force gmsh to use system C++ runtime (Ubuntu) rather than the ones overwritten in matlab env
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

# Optional: clear MATLAB’s library path if called from MATLAB
unset LD_LIBRARY_PATH

#/opt/homebrew/bin/
gmsh ${mesh_path}temp.geo -2 -format msh2 -o ${mesh_path}temp.msh 
