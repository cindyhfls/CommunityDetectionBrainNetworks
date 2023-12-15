% Script for running infomap and making models of the functional networks
% in some group of data
%% Initialization
clear;clc;close all;
% path_to_code = '/data/wheelock/data1/people/Muriah/code/NLA/Codes/FinalCode/BrBx-HSB_110418';
path_to_code = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup' % state the directory for this code
cd(path_to_code);
addpath(genpath(path_to_code));
infomappath = fullfile(path_to_code,'ExternalFunctions/infomap');
%% Load parameters

params = get_params_mat_eLABE_Y2_N92_healthyterm(infomappath);
% [params] = get_params_mat_WashU120(infomappath);
%% Use Evan's wrapper?

% Run_Infomap(params.zmat, params.dmat, params.xdist, params.lo:params.step:params.hi, params.binary,  params.outputdir, 2,[],500)

%% Wrapper to Infomap

stats=GraphCluster_HSB(params.zmat,params); % this can take 10-20 minutes depending on your data/sever specs
%% Save the results in a stats structure

stats.params.dmat = [];% let's not save the dmat as they might be hugeee for vertex-wise
stats.MuMat = [];% let's not save the matrix as they might be hugeee for vertex-wise

cd(params.outputdir)
save(fullfile(params.outputdir,params.IMap_fn),'stats'); % save infomap output to matrix


%% Calculate some stats
% stats.metrics = Matrix_metrics_HSB(stats.clusters,stats.MuMat,stats.rth,stats.params.binary);
% stats.metrics = Matrix_metrics_HSB_mod(stats);

%% Display something about results

if stats.params.fig
    stats.metrics.kdenth = stats.kdenth;stats.metrics.rth = stats.rth;
    Plot_InfoMap_Metrics_HSB(stats.metrics,stats.clusters); % Not quite clearly labeled   
    stats.metrics =rmfield(stats.metrics,{'kdenth','rth'});
end

%% Back up this file

fn = mfilename('fullpath');
mfilebackup(fn,params.IMap_fn);


