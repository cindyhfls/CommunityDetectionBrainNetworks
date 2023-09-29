function plot_network_assignment_parcel_key(Parcels, key,cmap,nets,removenone)
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
    cmap = linspecer(nNet);
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

% figure;
% if exist('nets','var')&& ~isempty(nets)   
%     ax = subplot(2,2,1);
%     params.fig_handle = ax;
%     params.view='lat';        % also, 'post','lat','med'
%     PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
%     title(sprintf('K = %i',nNet),'FontSize',15)
%     ax = subplot(2,2,2);
%     params.fig_handle = ax;
%     params.view='med';        % also, 'post','lat','med'
%     PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
%     currnets = nets(setdiff(unique(key),[0;none_idx]));
%     currcmap = cmap(setdiff(unique(key),[0;none_idx]),:);
%     
%     subplot(2,2,[3,4]);
%     h = gscatter(ones(1,nNet),ones(1,nNet),currnets,currcmap,'s',50);
%     for i = 1:nNet
%         set(h(i),'Color','k','MarkerFaceColor',currcmap(i,:));
%     end
%     legend(currnets,'interpreter','none','FontSize',9,'location','best','Orientation','horizontal','NumColumns',2);
%     xlim([10,11]);
%     axis('off')
% else
    ax = subplot(2,1,1);
    params.fig_handle = ax;
    params.view='lat';        % also, 'post','lat','med'
    PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
    title(sprintf('K = %i',nNet),'FontSize',15)
    ax = subplot(2,1,2);
    params.fig_handle = ax;
    params.view='med';        % also, 'post','lat','med'
    PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
% end
end