function Explore_parcel_kden_HSB(ROIclust,Cmap,Parcels,kden,fn)
%
% This function generates an interactive display of ROI Sortings for
% various kden, scanning through the kden chosen.
global keyPressed
fprintf('Press colormaps to switch thresholds\n, and press q for quit and s for saving figures when the mouse press is outside the colormaps\n')
if ~exist('fn','var')||isempty(fn)
    fn = 'tmp';
end
%% Set up paramters
[Nroi,Nkden]=size(ROIclust);
if any(ROIclust(:)==0)
    Cmap=cat(1,Cmap,[0.5,0.5,0.5]);
    ROIclust(ROIclust==0)=size(Cmap,1);
end
params.Cmap.P=Cmap;%IM.cMap;jet(nNet)
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'

button=0;
SortRowVec=[1:Nkden,1:Nkden,1:Nkden];
j=1;        % initial display


%% Display figures and write movie
f=figure('Position',[310,530,1550,570]);
while button~=2
 
%% Assign colors
load('MNI_coord_meshes_32k.mat');
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr
[Parcel_Nets.CtxL,Parcel_Nets.CtxR] = deal(NaN(size(Parcels.CtxL)));
key =ROIclust(:,j);
for ii = 1:size(key,1)
    Parcel_Nets.CtxL(Parcels.CtxL==ii,1) = key(ii);
    Parcel_Nets.CtxR(Parcels.CtxR==ii,1)= key(ii);
end
Anat.CtxL.data=Parcel_Nets.CtxL;
Anat.CtxR.data=Parcel_Nets.CtxR;
%% Plot
subplot(2,4,[1,5],'Position',[.025,0.025,.28,.9])
params.view= 'dorsal';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
title(['kden= ',num2str(kden(j),'%0.3f')],'Color','k')

subplot(2,4,[2:3],'Position',[.305,0.505,.44,.45])
params.view='lat';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

subplot(2,4,[6:7],'Position',[.305,0.005,.44,.45])
params.view='med';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

subplot(2,4,[4,8],'Position',[.775,0.01,.20,.95])
hold off
imagesc(sortrows(ROIclust,j));colormap(Cmap);hold on;axis off
plot([j-.5,j-0.5],[0,size(ROIclust,1)],'k');
plot([j+.5,j+0.5],[0,size(ROIclust,1)],'k');
set(gcf,'Color','w')
text(j-1,-2,'*','FontSize',20)    

% To better appearance, try organizing data a couple of rounds
foo=ROIclust;
SRV=circshift(SortRowVec,[(j-1),0]);
for k=1:(3*Nkden)
    foo=sortrows(foo,SRV(k));
end

foo(foo==0)=NaN;
image(sortrows(foo,j));colormap(Cmap);hold on;axis off
plot([j-.5,j-0.5],[0,size(ROIclust,1)],'k');
plot([j+.5,j+0.5],[0,size(ROIclust,1)],'k');
set(gcf,'Color','w')
text(j-0,-2,'*','FontSize',20)

%% get user input to change to a different threshold
[x,y,button]=ginput(1);

tmp=round(x);
if tmp<0 || tmp>size(ROIclust,2) 
    keyPressed = [];
    while isempty(keyPressed) % pressed outside the possible limits
        set(f, 'KeyPressFcn', @keyPressFunction);      
        pause(0.1)
    end
    if strcmp(keyPressed, 's')
        % Define your directory and filename
        fullPath = [fn,num2str(kden(j)),'.png'];
        % Save the figure
        saveas(gcf, fullPath);
        disp(['Figure saved as ' fullPath]);
    elseif strcmp(keyPressed, 'q')
        return
    end
else
    j = tmp;
end

end
end

function keyPressFunction(~,event)
global keyPressed
     keyPressed = event.Key;
end