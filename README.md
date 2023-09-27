# Infomap_MATLAB_wrapper for community detection in neuroimaging data

This package is currently under development and testing.

It takes a sparsely thresholded Nroi x Nroi structural/functional connectivity matrix and find the community structure using infomap community detection (mapequation.org).
Then it calculates the stability of solutions across algorithm initation, and across different resolution scale (by varying the threshold to the graph).

1. Code from other toolbox (External Functions):
   - infomap (mapequation.org) version 0.x 
   - BCT toolbox (2019_03_03_BCT version) download from https://sites.google.com/site/bctnet/ 
   - Faskowitz2018wsbmLifeSpan https://github.com/faskowit/Faskowitz2018wsbmLifeSpan/tree/master
   - NetworkCommunityToolbox
2. Currently support parcel-based/ROI-based input (Nroi x Nroi).
3. Full support for vertex-based and voxel-based input in progress.

Major updates:
2023/09/03: bringing distance exclusion after the disirable threshold and allow thresholding from backbone MST as a potential option. 
