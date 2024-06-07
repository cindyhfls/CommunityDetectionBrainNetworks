function View_Single_Assignment_Cortex(parceldata,cMap,parcelname)
Nparcels = length(setdiff(unique(parceldata),0))

if ~exist('cMap','var')
    cMap = jet(Nparcels);
end
if ~exist('parcelname','var')
    parcelname = '';
end


Linds=with_without_mw_conversion('Lindfull');
Rinds=with_without_mw_conversion('Rindfull');

n_verts_per_hem = 32492;
Parcels.CtxL = zeros(n_verts_per_hem,1);% total vertices
Parcels.CtxR = zeros(n_verts_per_hem,1);

Parcels.CtxL(Linds) = parceldata(1:length(Linds));
Parcels.CtxR(Rinds-n_verts_per_hem) = parceldata(length(Linds)+1:end);

%% View all parcels on MNI
load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr


% (1) Add parcel nodes to Cortices
% Parcels.CtxL(Parcels.CtxL==0) = Nparcels+1; % set the last one to black to draw borders (if available) and medial wall
% Parcels.CtxR(Parcels.CtxR==0) = Nparcels+1;
Anat.CtxL.data=Parcels.CtxL;
Anat.CtxR.data=Parcels.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=cMap;
% params.Cmap.P(Nparcels+1,:) = [0 0 0];
params.TC=1;
params.ctx='inf';           % 'std','inf','vinf'
figure;
t = tiledlayout(2, 1, "TileSpacing", "tight");
ax = nexttile;
% ax = subplot(2,1,1);
params.fig_handle = ax;
params.view= 'lat';       % 'dorsal','post','lat','med'
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
title(parcelname,'interpreter','none','color','k')
ax = nexttile;
% ax = subplot(2,1,2);
params.fig_handle = ax;
params.view ='med';
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

set(gcf,'color','w','InvertHardCopy','off');
end