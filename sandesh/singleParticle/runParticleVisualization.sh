#!/bin/bash

# Source OpenFOAM environment
# source /usr/lib/openfoam/openfoam2312/etc/bashrc  

# Build the utility
echo "Building particle visualization utility..."
wmake

# Generate mesh
echo "Generating mesh..."
blockMesh

# Run the particle visualization
echo "Running particle visualization..."
./particleVisualization

echo "Done! You can now open ParaView and load the case to visualize the particle." 