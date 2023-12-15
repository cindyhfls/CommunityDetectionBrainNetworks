function [params] = get_params_mat_WashU120(infomappath)
    %%
    params.numworkers = 0; % set to 0 to not use parallel
    if params.numworkers>0
        params.version = '0.15.7'; % Powers et al. 2011 Neuron version somehow the parallel acts weird in the 0.18.9 version but this saves a different clu file (only assignments no codelength)
    else
        params.version ='0.18.9'; % Eggebrecht et al. 2017 Cerebral Cortex
    end
    params.infomappath = fullfile(infomappath,sprintf('Infomap-%s',params.version),'Infomap');
    %% Basic information about my dataset
    outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap';

    params.format = 'mat'; % 'mat' or 'cifti' 
    params.parcel_name = 'Gordon';
    params.zmatfile = '/data/wheelock/data1/datasets/WashU120/pconns/washu120_parcellation_Gordon_20231101.mat'; % assumes the matrix is organized as nroi x nroi x nsubjects

    params.subjidx  = true(120,1);
    
    params.datasetname = 'WashU120' % folder name to store outputs
       
    timestr=datestr(datetime('now'),'yymmdd');
    
    params.outputdir = fullfile(outputdir,'parcel-wise',params.datasetname,params.parcel_name, timestr);
    if exist(params.outputdir,'dir') && (length(dir(params.outputdir))>2)
        params.outputdir = create_unique_directory(params.outputdir);
    end
    %% Load data    
    params.zmat = smartload(params.zmatfile);
    params.zmat = mean(params.zmat(:,:,params.subjidx),3);
    
    params.dmatfile = ['Parcels_',params.parcel_name,'.mat'];
    load(params.dmatfile,'parcels_dmat','ROIxyz');
    params.roi = ROIxyz;   % Coordinates for ROIs, used with exclusion distance
    params.dmat = parcels_dmat;% use geodesic distance
   %% Parameters that can vary
    params.repeats_consensus = 1; % whether or not to get all repeats
    params.binary=0;     % whether or not to binarize the matrix. Default=0;
    params.type = 'kden'; % choose between 'kden','r' and 'mst' for 'density threshold','raw correlation threshold','maximum spanning tree threshold'
    if strcmp(params.type,'mst')
        N = length(zmat);
        params.lo = 2/N; % use the MST density as minimum (N-1)/(N*(N-1)/2)
    else
        params.lo=0.010;      % Edge density minimum, typically 1% for ROIs
    end
    params.step=0.001;   % Edge density step, typically 0.001
    params.hi=0.20;      % Edge density maximum, typically 0.1
    params.xdist=20;     % Exclusion distance to minimize PSF shared variance
    params.fig = 0; % plot some figures
    if strcmp(params.format,'mat')
        params.repeats = 500; % default parameter assuming infomap convergence at n repeats, do NOT change unless you know what it is % Power et al. 2011 used 1000 but I think that's too many
    elseif strcmp(params.format,'cifti')
        params.repeats = 100; % default parameter assuming infomap convergence at n repeats, do NOT change unless you know what it is
    end   
    %% Make output directories
    if ~exist(params.outputdir,'dir')
        mkdir(params.outputdir);
    end
    params.IMap_fn=sprintf('Infomap_%s_low%1.3f_step%1.3f_high%1.3f_xdist%i.mat',params.datasetname,params.lo,params.step,params.hi,params.xdist); % name your output so you know what it is later
end