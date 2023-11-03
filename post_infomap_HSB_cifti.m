%% If loading from txt
clear
% params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/WashU120/20230502'

params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/eLABE_Y2/20230509'
stats.kdenth = importdata(fullfile(params.outputdir,'thresholds.txt'));

if ~isfield(stats,'removalidx')
    stats.removalidx = false(59412,1);
    if exist(fullfile(params.outputdir,'NaN_idx.mat'),'file')
        load(fullfile(params.outputdir,'NaN_idx.mat'),'idx')
        stats.removalidx(idx) = true;
    end
end

if sum(stats.removalidx)>0
    stats.clusters = zeros(length(stats.removalidx),size(stats.kdenth,2));
    stats.clusters(~stats.removalidx,:) = importdata(fullfile(params.outputdir,'rawassn.txt'));
end

%%
minsize = 400;
nameoption = 3;% 1: automatic, 3: using template
% stats.SortClus =OrgClustMat_HSB(stats.clusters,minsize,0); % last argument = 1 for sorting from higher threshold to lower threshold
stats.SortClus = remove_singleton(stats.clusters,minsize);
stats.SortClus = postprocess_ordinal_multilayer(stats.SortClus);
stats.SortClus = rename_multiscale(stats.SortClus);

templatepath  = 'Gordon2017_17Networks.dlabel.nii';%'Tu_eLABE_Y2_22Networks.nii'%'Gordon2017_17Networks.dlabel.nii';
ParcelCommunities =cifti_read(templatepath); % still use the Gordon colors for an unknown parcel?
[CW,GenOrder] =assign_Infomap_networks_by_template_cifti(stats.SortClus,ParcelCommunities,0.1,'dice');%'dice'
% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=stats.SortClus;foo(foo==0)=NaN;
stats.SortClusRO=stats.SortClus;
for j=1:length(GenOrder),stats.SortClusRO(foo==GenOrder(j))=j;end
load('Parcels_cortex_nomedialwall.mat')
%% (optional) Viewing and Manual edit of specific networks
for iNet =8:9
    Edit_NetworkColors(stats.SortClusRO,CWro,iNet,Parcels);
%     pause;
%     close all;
end
% customcolor = [];% fill in if you want to change it
% customcolor = distinguishable_colors(1,CWro.cMap); % set color to
% one
% CWro=Edit_NetworkColors(newstats.SortSortClus,CWro,iNet,Parcels,customcolor);
%%
Explore_parcel_kden_HSB(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'kden'));
%% Make video
Make_parcel_kden_Video(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'video'))
