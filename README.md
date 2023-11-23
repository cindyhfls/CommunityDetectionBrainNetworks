# Infomap_MATLAB_wrapper for community detection in neuroimaging data

## Introduction and motivation
This package is currently under development and testing. The initial commit is done by Dr. Adam Eggebrecht and is an upgrade from the [GraphTools] (https://www.nitrc.org/projects/graphtools/), with an attempt to make community detection of brain networks more principled and data-driven, especially in optimizing for the scale (number of communities).

It takes a sparsely thresholded Nroi x Nroi structural/functional connectivity matrix and find the community structure using infomap community detection (Rosvall & Bergstrom 2008 PNAS).
Then it calculates the stability of solutions across algorithm initation, and across different resolution scale (by varying the threshold to the graph). 

It is based on multiple lines of evidence which suggest a hierarchical organization in the brain while the original applications (e.g. Power 2011 Neuron) aimed to find a consensus across scales. I plan to determine the final communities in brain structural/functional connectivity data based on the assumption that at an optimal scale, 1) a given community detection algorithm is able to produce similar solution with each random initation, 2) the solution is similar across subsamples/bootstrapped samples of the input data, 3) the solution is similar across different subsamples/bootstraped samples of the input data, 4) clustering measures, such as silhouette index is relatively high, 5) the neighboring scales produce similar solutions and the transition is sharp across different hierarchies, 6) convergent evidence can be found from different community detection algorithm

The main goal for this repository is to explore the evidences above to find the optimal scales of human brain network organization and improves the visualization and user interaction/flexibility to run the infomap algorithm as well as alternative community detection algorithms including modularity maximization, spectral clustering, k-means clustering, weighted stochastic block model. In the future, I also intend to add detection of communities using multilayer networks to examine the linked community evolution across lifespan and/or states.

Additional goals include to correct for some bugs in the program and introduce the MST-based (Hagmann eta l. 2008 PLOS Biology) and local density based (Gordon et al. 2020 PNAS) threshold such that the graph is connected at the sparsest density level and produce less noisy solutions at those densities to enable the discovery of sub-networks beyond the commonly examined scale (e.g. in Gordon et al.,2020 PNAS, sub-networks within the default mode network were delineated). And to minimize the effect of user-defined parameters on the results by examining the effect of changing those parameters and provide suggestions of the default parameters, modularize the package to separate the community detection, evaluation and selection, and visualization.

## Working flow
- run_infomap_HSB.m loads data and saves the community detection output at different thresholded sparse graphs from the brain structural/functional connectivity matrix. This step is fully functional and saves the output solutions along with different diagnostic statistics.
- post_infomap_HSB*.m evaluates the different solutions to find the optimal set of consensus solutions and visualize it. Still developing and cleaning up.

## Notes:
1. Code from other toolbox (External Functions):
   - [infomap version 0.x](https://www.mapequation.org/)) 
   - [cifti-matlab](https://github.com/Washington-University/cifti-matlab)
   - [BCT toolbox 2019_03_03_BCT version](https://sites.google.com/site/bctnet/) 
   - [Faskowitz2018wsbmLifeSpan](https://github.com/faskowit/Faskowitz2018wsbmLifeSpan/tree/master)
   - [NetworkCommunityToolbox](http://commdetect.weebly.com/)
2. Currently support parcel-based/ROI-based input (Nroi x Nroi).
3. Full support for vertex-based and voxel-based input in progress.
