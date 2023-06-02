% Script for running infomap and making models of the functional networks
% in some group of data
clear;clc;close all;
path_to_code = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup' % state the directory for this code
cd(path_to_code);
addpath(genpath(path_to_code));

%% set your output directory to save files to
outputdir = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup/Results/';
params.format = 'mat'; % 'mat' or 'cifti' % N.B.: has not tested cifti yet
zmatfile = './ExampleData/120_allsubs_corr_Gordon.mat'; % file path for FC
datasetname = 'washu120'; % user note of the data
parcel_name = 'Gordon';% leave empty if params.format = cifti
outputdir = fullfile(outputdir,datasetname,datestr(datetime('now'),'yymmdd'));

%% Set parameters for use in Infomap
params.binary=0;     % Whether or not to Infomap with weights. Default=0;
params.lo=0.001;      % Edge density minimum, typically 1% for ROIs
params.step=0.001;   % Edge density step, typically 0.001
params.hi=0.1;      % Edge density maximum, typically 0.1
params.xdist=20;     % Exclusion distance to minimize PSF shared variance
params.fig = 0; % plot some figures
if strcmp(params.format,'mat')
    params.repeats = 1000; % default parameter assuming infomap convergence at n repeats, don't change unless you know what it is
    params.killTH = 5; % default parameter for the minimum number of nodes in the final consensus to be considered a network, this can be changed but the default is usually fine for group average data
elseif strcmp(params.format,'cifti')
    params.repeats = 100; % default parameter assuming infomap convergence at n repeats, don't change unless you know what it is
    params.killTH = 400; % default parameter for the minimum number of nodes in the final consensus to be considered a network, this can be changed but the default is usually fine for group average data
end
tmp = dir(zmatfile);
params.zmatfile = fullfile(tmp.folder,tmp.name);
timestr=datestr(datetime('now'),'yymmdd');
params.IMap_fn=sprintf('Infomap_%s_low%1.3f_step%1.3f_high%1.3f_xdist%i.mat',datasetname,params.lo,params.step,params.hi,params.xdist); % name your output so you know what it is later
params.outputdir = outputdir;
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
        params.roi = ROIxyz;   % Coordinates for ROIs, used with exclusion distance
        params.dmat = parcels_dmat;% use geodesic distance
end

%% Wrapper to Infomap % %

stats=GraphCluster_HSB(zmat,params); % this can take 10-20 minutes depending on your data/sever specs
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end

%% Display something about results
if stats.params.fig
    Plot_InfoMap_Metrics_HSB(stats); % Not quite clearly labeled   
end

%% Run Consensus Procedure to reduce the number of networks
% uses normalized mutual information across thresholds specified in params

[Cons,stats]=Org_Cons_Org_IMap_Matrix_HSB(stats);

%% Save the results in a stats structure

stats.params.dmat = [];% let's not save the dmat as they might be hugeee for vertex-wise
stats.MuMat = [];% let's not save the matrix as they might be hugeee for vertex-wise

cd(outputdir)
save(params.IMap_fn,'stats','Cons'); % save infomap output to matrix
% print(gcf,['./Figures/',params.IMap_fn,'_stats'],'-dtiff')

close all;



%% Back up this file
fn = mfilename('fullpath');
mfilebackup(fn,params.IMap_fn);


%% Supporting functions
% FileNameAndLocation=[mfilename('fullpath')];
function mfilebackup(FileNameAndLocation,identifier)
if ~exist('identifier','var')
    identifier = '';
end
timestr=datestr(datetime('now'),'yyyymmdd HH:MM:SS');
[filepath,name] = fileparts(FileNameAndLocation);
if ~exist(fullfile(filepath,'Backup'),'dir')
    mkdir(fullfile(filepath,'Backup'));
end
newbackup=fullfile(filepath,'Backup',sprintf('%s_%s_%s.m',name,timestr,identifier));
currentfile=strcat(fullfile(filepath,[name, '.m']));
copyfile(currentfile,newbackup);
end

