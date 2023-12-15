function CW = Edit_NetworkColors(Clust,CW,iNet,Parcels,customcolor)
if exist('customcolor','var')&&~isempty(customcolor)
    assert(all(size(customcolor)==[1,3]))
    CW.cMap(iNet,:) = customcolor;
end
load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr
%% Select largest representation of a module label
G1=setdiff(unique(Clust(:)),0);
for j=1:length(G1)
    tmp =sum(Clust==G1(j,1));tmp(tmp==0) = NaN;
    [G1(j,2)]=mode(tmp); % Adam used max when making manual judgement of network names
    G1(j,3) = find(tmp==G1(j,2),1);
end

%% View only parcels within a specific network
Parcel_Nets.CtxL = double(any(Parcels.CtxL==find(Clust(:,G1(iNet,3)) == G1(iNet,1))',2))*iNet; 
Parcel_Nets.CtxR = double(any(Parcels.CtxR==find(Clust(:,G1(iNet,3)) == G1(iNet,1))',2))*iNet;

Anat.CtxL.data=Parcel_Nets.CtxL; % plot desired networks
Anat.CtxR.data=Parcel_Nets.CtxR;

keep = false(length(CW.Nets),1);
keep(iNet) = true;
cMap=CW.cMap;
cMap(~keep,:)=0.5; % color to gray for other networks

figure('position',[100 100 800 800]);
params.Cmap.P=cMap;
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'
params.view='lat';        % also, 'post','lat','med','dorsal'
params.fig_handle = subplot(2,1,1);
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
title(['N',num2str(iNet),': ',CW.Nets{iNet}],'interpreter','None')
params.fig_handle = subplot(2,1,2);
params.view = 'med';
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

end