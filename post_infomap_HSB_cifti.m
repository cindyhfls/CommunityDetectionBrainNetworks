%% If loading from txt
clear
% params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/WashU120/20230502'

% params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/eLABE_Y2/20230509'
params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/eLABE_Y2_N92/20231102'
stats.kdenth = importdata(fullfile(params.outputdir,'thresholds.txt'));
load(fullfile(params.outputdir,'params.mat'));

if ~isfield(params,'removalidx')
    params.removalidx = false(59412,1);
    if exist(fullfile(params.outputdir,'NaN_idx.mat'),'file')
        load(fullfile(params.outputdir,'NaN_idx.mat'),'idx')
        params.removalidx(idx) = true;
    end
end

if sum(params.removalidx)>0
    stats.clusters = zeros(length(params.removalidx),size(stats.kdenth,2));
    stats.clusters(~params.removalidx,:) = importdata(fullfile(params.outputdir,'rawassn.txt'));
end

%%
minsize = 400;
nameoption = 3;% 1: automatic, 3: using template
stats.SortClus   = postprocess_ordinal_multilayer(stats.clusters);  % sort across level to have some consistency across columns
stats.SortClus = remove_singleton(stats.SortClus ,minsize);
templatepath  = 'Gordon2017_17Networks.dlabel.nii';%'Tu_eLABE_Y2_22Networks.nii'%'Gordon2017_17Networks.dlabel.nii';
ParcelCommunities =cifti_read(templatepath); % still use the Gordon colors for an unknown parcel?
[CW,GenOrder] = assign_networks_by_template_cifti(stats.SortClus,ParcelCommunities,0.1,'dice')

% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=stats.SortClus;foo(foo==0)=NaN;
stats.SortClusRO=stats.SortClus;
for j=1:length(GenOrder),stats.SortClusRO(foo==GenOrder(j))=j;end
load('Parcels_cortex_nomedialwall.mat')
%% (optional) Viewing and Manual edit of specific networks
for iNet =9:11
    Edit_NetworkColors(stats.SortClusRO,CWro,iNet,Parcels);
%     pause;
%     close all;
end
% customcolor = [];% fill in if you want to change it
% customcolor = distinguishable_colors(1,CWro.cMap); % set color to
% one
% CWro=Edit_NetworkColors(newstats.SortSortClus,CWro,iNet,Parcels,customcolor);
%% Save solution
cd(params.outputdir)
save(['solution_minsize_',num2str(minsize),'.mat'],'CWro','stats');

%%
Explore_parcel_levels_HSB(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'kden'),40);

%% Make video
Make_parcel_kden_Video(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'video'))
