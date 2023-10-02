# Infomap_MATLAB_wrapper for community detection in neuroimaging data

This package is currently under development and testing. The initial commit is done by Dr. Adam Eggebrecht and is an upgrade from the GraphTools https://www.nitrc.org/projects/graphtools/, with an attempt to make community detection of brain networks more principled and data-driven, especially in optimizing for the scale (number of communities).

It takes a sparsely thresholded Nroi x Nroi structural/functional connectivity matrix and find the community structure using infomap community detection (mapequation.org).
Then it calculates the stability of solutions across algorithm initation, and across different resolution scale (by varying the threshold to the graph). 

It is based on multiple lines of evidence which suggest a hierarchical organization in the brain while the original applications (e.g. Power 2011 Neuron) aimed to find a consensus across scales. I plan to determine the final communities in brain structural/functional connectivity data based on the assumption that at an optimal scale, 1) a given community detection algorithm is able to produce similar solution with each random initation, 2) the solution is similar across subsamples/bootstrapped samples of the input data, 3) at the 
and potentially across different subsamples/bootstraped samples of the input data, 4) clustering measures, such as silhouette index is relatively high, 5) the neighboring scales produce similar solutions and the transition is sharp across different hierarchies, 6) convergent evidence can be found from different community detection algorithm

The goal for this repository is to explore the evidences above to find the optimal scales of human brain network organization and improves the visualization and user interaction/flexibility to run the infomap algorithm as well as alternative community detection algorithms including modularity maximization, spectral clustering, k-means clustering, weighted stochastic block model. In the future, I also intend to add detection of communities using multilayer networks to examine the linked community evolution across lifespan and/or states.

Notes:
1. Code from other toolbox (External Functions):
   - infomap (mapequation.org) version 0.x 
   - BCT toolbox (2019_03_03_BCT version) download from https://sites.google.com/site/bctnet/ 
   - Faskowitz2018wsbmLifeSpan https://github.com/faskowit/Faskowitz2018wsbmLifeSpan/tree/master
   - NetworkCommunityToolbox
2. Currently support parcel-based/ROI-based input (Nroi x Nroi).
3. Full support for vertex-based and voxel-based input in progress.

Major updates:
2023/09/03: bringing distance exclusion after the disirable threshold and allow thresholding from backbone MST as a potential option. 
