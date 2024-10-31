# Community Detection in Brain Networks

Functional brain networks can be found with unsupervised, data-driven community detection/clustering methods using functional connectivity data. Sometimes researchers from network science/statistics backgrounds fail to recognize the biological and functional relevance of the brain networks, and obtain solutions purely based on algorithms and lack biological interpretation. On the other hand, researchers from biological/medical backgrounds would run a selected algorithm assuming that it will just produce the desirable outcome without recognizing the stochasticity in the algorithm, the stability of their solution, and the selection of the number of clusters or consensus from a group of solutions. As a result, different network divisions were obtained from the same open dataset (e.g. Tooley et al. 2022 Neuroimage and Marek et al. 2019 Neuroimage on the Adolescent Brain Cognitive Development).

Here, my goal is to provide an integrated toolbox that not only considers different community detection algorithms but also focuses on "what next" after one gets the solutions. In particular, I want to guide people through the recent advances in network science in solving the degeneracy problem of community detection, identifying the most likely number of clusters (or a hierarchy of them), and the intuition behind different community detection algorithms or consensus methods. In this way, researchers and audiences can better interpret their results with respect to the broader literature.

One additional confound is the nomenclature of the networks (Uddin et al. 2023 Network Neuroscience; Uddin et al. 2019 Brain topography) and the elusive roles of putative "clusters" identified in resting-state fMRI connectivity without validation with task fMRI activation using functional localizers, and therefore do not have a meaningful interpretation of their functional role. Networks may serve diverse roles and cognitive tasks usually depend on the involvement of different cognitive abilities. Designing an experiment that separates a single non-sensory function is non-trivial. I want to use the existing meta-analysis tools (https://neurosynth.org/ and https://nimare.readthedocs.io/en/stable/) as well as their overlap in topography and connectivity profiles to existing networks (e.g. Network Correspondence Toolbox https://www.biorxiv.org/content/10.1101/2024.06.17.599426v1) to help make sense of and naming the "clusters" identified using connectivity data only.

N.B. In addition to the purely data-driven approaches, there are also recent developments in obtaining individualized functional networks based on a group prior (e.g. Cui et al. 2020 Neuron, Kong et al. 2019 Cerebral Cortex; Hacker et al. 2013 Neuroimage; Gordon et al. 2017 Cerebral Cortex). I might briefly discuss it at the end. 

This is developed with human neuroimaging data (fMRI) in mind but has the potential to generalize to other modalities of neuroimaging and animal studies.

## Original version and planned updates
This toolbox is currently under development and testing. The initial commit was the version written by Dr. Adam Eggebrecht (Eggebrecht et al. 2017 Cerebral Cortex) as an upgrade from the [GraphTools] (https://www.nitrc.org/projects/graphtools/). 

It takes a sparsely thresholded Nroi x Nroi structural/functional connectivity matrix and find the community structure using infomap community detection (Rosvall & Bergstrom 2008 PNAS).
Then it calculates the stability of solutions across algorithm initiation, and across different resolution scales (by varying the threshold to the graph). This approach was initially introduced in Power et al. 2011 Neuron. 

Of the multiple ways to determine the number of clusters, I plan to use a combination of validations e.g. 1) a given community detection algorithm can produce similar solutions with each random initiation, 2) the solution is similar across subsamples/bootstrapped samples of the input data, 3) the solution is similar across different subsamples/bootstrapped samples of the input data, 4) clustering measures, such as silhouette index is relatively high, 5) the neighboring scales produce similar solutions and the transition is sharp across different hierarchies, 6) convergent evidence can be found from different community detection algorithms

I would explore the evidence above to find the optimal scales of human brain network organization and improve the visualization and user interface/flexibility to run the infomap algorithm as well as alternative community detection algorithms including modularity maximization, spectral clustering, k-means clustering, weighted stochastic block model. In the future, I also intend to add the detection of communities using multilayer networks to examine the linked community evolution across lifespans and/or states.

Additional goals include correcting for some bugs in the original version of the Infomap community detection on sparsified connectivity matrix (Eggebrecht et al. 2017 Cerebral Cortex; Power et al. 2011 Neuron) and introducing the MST-based (Hagmann et al. 2008 PLOS Biology) and local density-based (Gordon et al. 2020 PNAS) threshold such that the graph is connected at the sparsest density level and producing less noisy solutions at those densities to enable the discovery of sub-networks beyond the commonly examined scale (e.g. in Gordon et al.,2020 PNAS, sub-networks within the default mode network were delineated). To minimize the effect of user-defined parameters on the results by examining the effect of changing those parameters and providing suggestions of the default parameters, set the random seed to obtain reproducible results and modularize the package to separate the community detection, evaluation and selection of the number of clusters, and visualization (both abstract graph-like visualization and network topography on the brain).

I also considered alternative null models that are more suitable for correlation matrix (Random-Matrix-Based filtering, MacMahon & Garlaschelli 2015 Physical Review), fitting distance-based thresholding instead of the current hard cut-off at 20 mm (Esfahlani et al. 2020 Neuroimage). For obtaining consensus from a group of solutions at different hierarchy, I consider adding the hierarchical consensus (Jeub et al. 2018 Scientific Reports).

## Tutorial:
post_comm_assignment_ordinal.mlx Interactive tutorial for post-community detection for different resolution scales or different temporal points
post_comm_assignment_categorical.mlx Interactive tutorial for post-community detection for different repeats/bootstraps/subjects/methods

## Notes:
1. Code from external toolboxes (External Functions):
   - [infomap (version 0.15.7 and 0.18.9)](https://www.mapequation.org/)) 
   - [cifti-matlab-master](https://github.com/Washington-University/cifti-matlab)
   - [BCT (2019_03_03 version)](https://sites.google.com/site/bctnet/) 
   - [Faskowitz2018wsbmLifeSpan](https://github.com/faskowit/Faskowitz2018wsbmLifeSpan/tree/master)
   - [NetworkCommunityToolbox](http://commdetect.weebly.com/)
   - [communityalg](https://github.com/CarloNicolini/communityalg/tree/master) (N.B.: downloaded some dependencies from here: https://github.com/CarloNicolini/algonet/tree/master)
   - [HierarchicalConsensus](https://github.com/LJeub/HierarchicalConsensus)
   - [hierarchical-brain-networks](https://github.com/emergelab/hierarchical-brain-networks) 
   - [GenLouvain](https://github.com/GenLouvain/GenLouvain)
   - [multilouvain](https://github.com/CarloNicolini/multilouvain)[origional C++ code from Vincent Traag](https://github.com/vtraag/louvain-igraph) 
2. Currently support parcel-based/ROI-based input (Nroi x Nroi).
3. Full support for vertex-based and voxel-based input is in progress.
