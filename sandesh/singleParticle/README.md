# Single Particle Visualization in OpenFOAM

This case demonstrates how to visualize a single particle in OpenFOAM using ParaView.

## Files Structure
```
singleParticle/
├── 0/
│   ├── positions      # Particle position (0.5, 0.5, 0.5)
│   ├── velocities     # Particle velocity (0.1, 0, 0)
│   └── diameters      # Particle diameter (0.1)
├── system/
│   ├── blockMeshDict  # Mesh definition (1x1x1 cube)
│   └── controlDict    # Control parameters
├── particleVisualization.C  # Main utility code
├── Make/
│   ├── files          # Build configuration
│   └── options        # Compiler options
├── runParticleVisualization.sh  # Run script
└── README.md          # This file
```

## Step-by-Step Instructions

### 1. Build and Run
```bash
# Make the run script executable
chmod +x runParticleVisualization.sh

# Run the complete process
./runParticleVisualization.sh
```

### 2. Manual Steps (if needed)
```bash
# Source OpenFOAM environment
source /usr/lib/openfoam/openfoam2312/etc/bashrc

# Build the utility
wmake

# Generate mesh
blockMesh

# Run visualization
./particleVisualization
```

### 3. Visualize in ParaView
1. Open ParaView
2. File → Open → Navigate to `singleParticle` folder
3. Select the case file and click "OK"
4. In the pipeline browser, you should see:
   - The mesh (blockMesh)
   - Lagrangian clouds (particleCloud)
5. Click "Apply" to load the data
6. The particle should appear as a sphere at position (0.5, 0.5, 0.5)

## Particle Properties
- **Position**: (0.5, 0.5, 0.5) - center of the domain
- **Velocity**: (0.1, 0, 0) - moving in x-direction
- **Diameter**: 0.1 - 10% of domain size

## Next Steps
Once you can visualize this single particle, we can proceed to:
1. Create a DEM particle class
2. Implement particle motion
3. Add multiple particles
4. Implement contact detection and forces

## Troubleshooting
- If ParaView doesn't show the particle, check that the lagrangian clouds are enabled
- If the utility doesn't build, ensure OpenFOAM environment is sourced
- Check the console output for any error messages 