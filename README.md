# astro-viskit

a collection of tools for analyzing and visualizing astronomical data.

## adaptive-histogram

Make 2d/3d histograms with adaptive grid refinement.

Refinement is not implemented using quad-tree/oct-tree, suitable for large data sets.

Function evaluation supported. For example, you can easily find the mean l.o.s velocity of stars in the histogram of [lon, lat]

## easy-histogram

Make 2d/3d histograms with function evaluation. Faster than adaptive-histogram, if uniform grid resolution is desired.

## voronoi-plot

Visualize data organized on Voronoi meshes. Originally designed to visualize kinematic maps of IFU datasets, now modified for general purposes.

*Other new tools will be migrated soon."
