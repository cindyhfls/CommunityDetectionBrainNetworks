function params = get_params_cifti_eLABE_Y2()
% place-holder now, need to change a lot of things here
    %% Basic information about my dataset
    outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap';

    params.format = 'cifti'; % 'mat' or 'cifti' % N.B.: has not tested cifti yet=
    params.parcel_name = '';
    params.zmatfile = '/data/wheelock/data1/datasets/eLABE/pconns/eLABE_Y2_N113_atleast600frames_parcellation_Tu_342_20231016.mat'
    
    healthysubj = importdata('/data/wheelock/data1/people/Cindy/BCP/subject_lists/eLABE_Y2_N92_healthyterm.txt');
    allsubj = importdata('/data/wheelock/data1/people/Cindy/BCP/subject_lists/eLABE_Y2_N113_atleast600frames.txt');
    subjidx = contains(allsubj,healthysubj);
    
    params.datasetname = 'eLABE_Y2_N92_healthyterm'
    params.subjidx = subjidx;
    
    timestr=datestr(datetime('now'),'yymmdd');
    
    if strcmp(params.format,'mat')
        outputdir = fullfile(outputdir,'parcel-wise',params.datasetname,params.parcel_name, timestr);
        params.parcel_name = parcel_name;
    elseif strcmp(params.format,'cifti')
        outputdir = fullfile(outputdir,'vertex-wise',params.datasetname, timestr);
    else
        error('Currently only supports mat and cifti format');
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    params.IMap_fn=sprintf('Infomap_%s_low%1.3f_step%1.3f_high%1.3f_xdist%i.mat',datasetname,params.lo,params.step,params.hi,params.xdist); % name your output so you know what it is later
    params.outputdir = outputdir;
   %% parameters that can vary
    params.repeats_consensus = 0;
    params.binary=0;     % Whether or not to Infomap with weights. Default=0;
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
    %% Load data
    switch params.format
        case 'cifti'
            % Get the fsLR 32k geodesic distance
            load('/data/nil-bluearc/GMT/Evan/Atlases/32k_ConteAtlas_v2_distribute/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat') % let's copy this file over later
            params.dmat = distances(1:59412,1:59412);
            clear distances
            
            % Get the FC matrix
            zmat =cifti_read(zmatfile);
            zmat = zmat.cdata;
            clear data;
        case 'mat'
            load(['./Parcels/Parcels_',parcel_name,'.mat'],'parcels_dmat','ROIxyz');
            %         load('Parcels_LR_distances.mat'); parcels_dmat = parcel_distances;% overwrite with the distance matrix that Evan had
            %         warning('please remove the line above after testing!');
            zmat = smartload(zmatfile);
            zmat = mean(zmat,3);
            params.roi = ROIxyz;   % Coordinates for ROIs, used with exclusion distance
            params.dmat = parcels_dmat;% use geodesic distance
    end
end