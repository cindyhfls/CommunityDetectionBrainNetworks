% Script for running infomap and making models of the functional networks
% in some group of data
%% Initialization
clear;clc;close all;
path_to_code = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup' % state the directory for this code
cd(path_to_code);
addpath(genpath(path_to_code));
global infomappath
infomappath =fullfile(path_to_code,'/ExternalFunctions/infomap/Infomap')
%% set your output directory to save files to
outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap';
params.format = 'mat'; % 'mat' or 'cifti' % N.B.: has not tested cifti yet
zmatfile = './ExampleData/120_allsubs_corr_Gordon.mat'; % file path for FC
% zmatfile = '/data/wheelock/data1/datasets/BCP/December2020/pconns/BCP_Dec2020_N177_parcellation_eLABE_Y2_prelim_05062023_20230616.mat' %'/data/wheelock/data1/datasets/eLABE/pconns/eLABE_Y2_N113_atleast600frames_parcellation_Gordon_20230530.mat'
datasetname ='WashU120'% 'BCP_Dec_N177'%'eLABE_Y2_N113'; % user note of the data
parcel_name ='Gordon'% leave empty if params.format = cifti %should match the name in folder ./Parcels/

if strcmp(params.format,'mat')
    outputdir = fullfile(outputdir,'parcel-wise',datasetname,parcel_name,datestr(datetime('now'),'yymmdd'));
    params.parcel_name = parcel_name;
elseif strcmp(params.format,'cifti')
    outputdir = fullfile(outputdir,'vertex-wise',datasetname,datestr(datetime('now'),'yymmdd'));
end
%% Load data
switch params.format
    case 'cifti'
        % Get the fsLR 32k geodesic distance
        load('/data/nil-bluearc/GMT/Evan/Atlases/32k_ConteAtlas_v2_distribute/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat')
        params.dmat = distances(1:59412,1:59412);
        clear distances
        
        lROI = gifti('/data/wheelock/data1/parcellations/SurfaceFiles/Conte69.L.sphere.32k_fs_LR.coord.gii');
        rROI = gifti('/data/wheelock/data1/parcellations/SurfaceFiles/Conte69.R.sphere.32k_fs_LR.coord.gii');
        ROIxyz = [lROI.vertices;rROI.vertices];
        params.roi = ROIxyz([Util.with_without_mw_conversion('Lindfull');Util.with_without_mw_conversion('Rindfull')],:);
        clear ROIxyz lROI rROI;
        % Get the FC matrix
        zmat =ft_read_cifti_mod(zmatfile);
        zmat = zmat.data;
        clear data;
    case 'mat'
        load(['./Parcels/Parcels_',parcel_name,'.mat'],'parcels_dmat','ROIxyz');
        zmat = smartload(zmatfile);
        zmat = mean(zmat,3);
        params.roi = ROIxyz;   % Coordinates for ROIs, used with exclusion distance
        params.dmat = parcels_dmat;% use geodesic distance
end
%% Set parameters for use in Infomap
params.repeats_consensus = 1;
params.binary=0;     % Whether or not to Infomap with weights. Default=0;
params.type = 'mst'; % choose between 'kden','r' and 'mst' for 'density threshold','raw correlation threshold','maximum spanning tree threshold'
if strcmp(params.type,'mst')
    [MST] =backbone_wu_mod(zmat);
    params.lo = nnz(MST)/length(zmat)^2; % use the MST density as minimum
else
    params.lo=0.001;      % Edge density minimum, typically 1% for ROIs
end
params.step=0.0025;   % Edge density step, typically 0.001
params.hi=0.15;      % Edge density maximum, typically 0.1
params.xdist=0;     % Exclusion distance to minimize PSF shared variance
params.fig = 0; % plot some figures
if strcmp(params.format,'mat')
    params.repeats = 500; % default parameter assuming infomap convergence at n repeats, do NOT change unless you know what it is % Power et al. 2011 used 1000 but I think that's too many
%     params.killTH = 5; % default parameter for the minimum number of nodes in the final consensus to be considered a network, this can be changed but the default is usually fine for group average data
elseif strcmp(params.format,'cifti')
    params.repeats = 100; % default parameter assuming infomap convergence at n repeats, do NOT change unless you know what it is
%     params.killTH = 400; % default parameter for the minimum number of nodes in the final consensus to be considered a network, this can be changed but the default is usually fine for group average data
end
tmp = dir(zmatfile);
params.zmatfile = fullfile(tmp.folder,tmp.name);
timestr=datestr(datetime('now'),'yymmdd');
params.IMap_fn=sprintf('Infomap_%s_low%1.3f_step%1.3f_high%1.3f_xdist%i.mat',datasetname,params.lo,params.step,params.hi,params.xdist); % name your output so you know what it is later
params.outputdir = outputdir;

%% Wrapper to Infomap % %

stats=GraphCluster_HSB(zmat,params); % this can take 10-20 minutes depending on your data/sever specs
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

%% Calculate some stats
% stats.metrics = Matrix_metrics_HSB(stats.clusters,stats.MuMat,stats.rth,stats.params.binary);
stats.metrics = Matrix_metrics_HSB_mod(stats);
%% Save the results in a stats structure

stats.params.dmat = [];% let's not save the dmat as they might be hugeee for vertex-wise
stats.MuMat = [];% let's not save the matrix as they might be hugeee for vertex-wise

cd(outputdir)
save(fullfile(outputdir,params.IMap_fn),'stats'); % save infomap output to matrix

close all;
%% Display something about results

if stats.params.fig
    stats.metrics.kdenth = stats.kdenth;stats.metrics.rth = stats.rth;
    Plot_InfoMap_Metrics_HSB(stats.metrics,stats.clusters); % Not quite clearly labeled   
    stats.metrics =rmfield(stats.metrics,{'kdenth','rth'});
end

%% Back up this file

fn = mfilename('fullpath');
mfilebackup(fn,params.IMap_fn);


