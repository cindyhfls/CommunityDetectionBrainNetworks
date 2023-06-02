function Explore_ROI_kden_HSB(ROIclust,Cmap,Anat,ROIxyz,kden)
%
% This function generates an interactive display of ROI Sortings for
% various kden, scanning through the kden chosen.

%% Set up paramters
[Nroi,Nkden]=size(ROIclust);
if any(ROIclust(:)==0)
Cmap=cat(1,Cmap,[0,0,0]);
ROIclust(ROIclust==0)=size(Cmap,1);
end

Params.Cmap=Cmap;
Params.Scale=length(Cmap);Params.OL=1;Params.PD=1;Params.TC=1;
Params.Th.P=0.001;Params.Th.N=-Params.Th.P;
ROI2.coord=ROIxyz;
ROI.radius=repmat(10,[Nroi,1]);
ROI.Nfaces=29;
Anat.ctx='inf';
Anat.alpha=1;
button=0;
SortRowVec=[1:Nkden,1:Nkden,1:Nkden];
j=1;        % initial display


%% Display figures and write movie
f=figure('Position',[310,530,1550,570]);
while button~=2
    
ROI2.color=Cmap(ROIclust(:,j),:);
if isfield(ROI2,'Network'),ROI2=rmfield(ROI2,'Network');end
ROI2.Network(:,1)=1:size(ROI2.coord,1);
ROI2.Network(:,2)=ones(size(ROI2.coord,1),1);

subplot(2,4,[1,5],'Position',[.025,0.025,.28,.9])
Anat.view='dorsal';Draw_ROIs_on_Cortex_HSB(Anat,ROI2);
title(['kden= ',num2str(kden(j),'%0.3f')],'Color','k')

subplot(2,4,[2:3],'Position',[.305,0.505,.44,.45])
Anat.view='lat';Anat.alpha=1;Draw_ROIs_on_Cortex_HSB(Anat,ROI2);

subplot(2,4,[6:7],'Position',[.305,0.005,.44,.45])
Anat.view='med';Anat.alpha=1;Draw_ROIs_on_Cortex_HSB(Anat,ROI2);

subplot(2,4,[4,8],'Position',[.775,0.01,.20,.95])
hold off

% To better appearance, try organizing data a couple of rounds
foo=ROIclust;
SRV=circshift(SortRowVec,[(j-1),0]);
for k=1:(3*Nkden)
    foo=sortrows(foo,SRV(k));
end

foo(foo==0)=NaN;
image(sortrows(foo,j));colormap(Params.Cmap);hold on;axis off
plot([j-.5,j-0.5],[0,size(ROIclust,1)],'k');
plot([j+.5,j+0.5],[0,size(ROIclust,1)],'k');
set(gcf,'Color','w')
text(j-0,-2,'*','FontSize',20)
[x,y,button]=ginput(1);
j=round(x);
end
