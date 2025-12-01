#!/usr/bin/bash

CANDIDATES=(
    "/opt/homebrew/bin/gmsh"
    "/usr/local/bin/gmsh"
    "/usr/bin/gmsh"
)

for p in "${CANDIDATES[@]}"; do
    if [ -x "$p" ]; then
        GMSH="$p"
        break
    fi
done

# Fallback to PATH search
if [ -z "$GMSH" ]; then
    GMSH=$(command -v gmsh)
fi

if [ -z "$GMSH" ]; then
    echo "Error: gmsh not found."
    exit 1
fi

echo "Using gmsh at: $GMSH"


mesh_path=$1

# Force gmsh to use system C++ runtime (Ubuntu) rather than the ones overwritten in matlab env
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

# Optional: clear MATLAB’s library path if called from MATLAB
unset LD_LIBRARY_PATH

$GMSH ${mesh_path}temp.geo -2 -format msh2 -o ${mesh_path}temp.msh 

#./mesh/msh2xdmf.py
