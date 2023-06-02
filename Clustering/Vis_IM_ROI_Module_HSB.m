function Vis_IM_ROI_Module_HSB(Clust,stats,Anat,k,Nroi)
%
% This function visualizes a single set of ROIs, grouped via infomap to aid
% in naming and coloring.
% Clust are the Nroi-x-Nkden clustering labels
% stats is the structure containing parameters for InfoMapping
% Anat contains info about the anatomy
% k is the network label of interest

%% Select largest representation of a module label
G1=unique(Clust(:));
for j=1:length(G1)
    [G1(j,2),G1(j,3)]=max(sum(Clust==G1(j,1)));
end

%% Set ROI info
ROI2.coord=stats.params.roi(Clust(:,G1(k,3))==G1(k,1),:);   % xyz
ROI2.color=repmat([1,0,0],size(ROI2.coord,1),1);            % color

ROI2.Network(:,1)=1:size(ROI2.coord,1);                     % Network
ROI2.Network(:,2)=ones(size(ROI2.coord,1),1);
if ~isfield(ROI2,'radius'), ROI2.radius=repmat(3,[Nroi,1]);end
ROI2.Nfaces=100;

%% Create Figure
figure('Position',[680,516,630,580],'visible','on');

subplot(3,6,[1,2,3,7,8,9])             % Dorsal solid view
Anat.view='dorsal'; 
Anat.alpha=1;
Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
title(['Col ',num2str(G1(k,3)),'; kden ',...
    num2str(stats.kdenth(G1(k,3)))],'Color','w')

subplot(3,6,[4,5,6,10,11,12]);          % Dorsal transparent view
Anat.alpha=.1;
Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
title(['Net ',num2str(G1(k,1))],'Color','w')

subplot(3,6,[13:14])                    % Left solid view
Anat.alpha=1;
Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
view([-90,0])

subplot(3,6,[15:16])                    % Right solid view
Anat.alpha=1;
Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
view([90,0])

subplot(3,6,[17:18]);                   % Lateral transparent view
Anat.alpha=.1;
Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
view([-90,0])





