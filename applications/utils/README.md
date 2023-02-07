# Pre- and post-processing utilities for dissolFoam

## surfRoughGen: 
OpenFOAM utility to generate a rough surface by Fourier synthesis.

usage: `surfRoughGen`

It has a dictionary `system/surfRoughGenDict` - the case `dissolFrac` uses `surfRoughGen`

## runMeshUpdateOnce: 
OpenFOAM utility for a single mesh update cycle.

usage: `runMeshUpdateOnce`

## orderBoundaries: 
User library adds support for rearranging boundaries after meshing to place moving boundary first.

