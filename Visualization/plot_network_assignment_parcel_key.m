function plot_network_assignment_parcel_key(Parcels, key,cmap,nets,removenone)
warning('For assignment in vertices it would be more efficient to use View_Single_Assignment_Cortex.m')
%% set defaults
if ~exist('removenone','var')||isempty(removenone)
   removenone =1; % don't count the network with the label 'None' as an actual network
end

if removenone && exist('nets','var') && ~isempty(nets)
    none_idx = find(string(nets)=='None');
    nNet = length(setdiff(unique(key),[0;none_idx]));
else
    nNet = length(setdiff(unique(key),[0]));
end
if ~exist('cmap','var')||isempty(cmap)
%     cmap = linspecer(nNet);
    cmap = distinguishable_colors(nNet);
end
%% find network assignments for each ROI
[Parcel_Nets.CtxL,Parcel_Nets.CtxR] = deal(NaN(size(Parcels.CtxL)));

for ii = 1:size(key,1)
    Parcel_Nets.CtxL(Parcels.CtxL==ii,1) = key(ii);
    Parcel_Nets.CtxR(Parcels.CtxR==ii,1)= key(ii);
end
%% Plot on inflated surface
load('MNI_coord_meshes_32k.mat');
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr
Anat.CtxL.data=Parcel_Nets.CtxL;
Anat.CtxR.data=Parcel_Nets.CtxR;
params.Cmap.P=cmap;%IM.cMap;jet(nNet)
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'

f = figure('position',[100 100 385 275]);
ax1 = subplot(2,1,1);
set(ax1,'Position',[0,0.5,1,0.5]);
params.lighting = 'None';
params.fig_handle = ax1;
params.view= 'lat';       % 'dorsal','post','lat','med'
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

ax2 = subplot(2,1,2);
set(ax2,'Position',[0,0,1,0.5]);
params.lighting = 'None';
params.fig_handle = ax2;
params.view ='med';
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);



end