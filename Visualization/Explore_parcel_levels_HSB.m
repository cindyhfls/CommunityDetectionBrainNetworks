function Explore_parcel_levels_HSB(ROIclust,Cmap,Parcels,levels,fn)
%
% This function generates an interactive display of ROI Sortings for
% various levels, scanning through the kden chosen.
global keyPressed
fprintf('Press on the block to switch between levels\n, and press q for quit and s for saving figures when the mouse press is outside the blocks\n')
if ~exist('fn','var')||isempty(fn)
    fn = 'tmp';
end
%% Set up paramters
[Nroi,Nlevels]=size(ROIclust);
if any(ROIclust(:)==0)
    Cmap=cat(1,Cmap,[0.5,0.5,0.5]);
    ROIclust(ROIclust==0)=size(Cmap,1);
end
params.Cmap.P=Cmap;%IM.cMap;jet(nNet)
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'

button=0;
SortRowVec=[1:Nlevels,1:Nlevels,1:Nlevels];
j=1;        % initial display


%% Display figures and write movie
f=figure('Position',[310,530,1550,570]);
while button~=2
    clf
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
subplot('Position',[.025,0.025,.28,.9])
params.view= 'dorsal';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
title(sprintf('%0.3f',levels(j)),'Color','k')

subplot('Position',[.305,0.505,.44,.45])
params.view='lat';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

subplot('Position',[.305,0.005,.44,.45])
params.view='med';
params.fig_handle = gca;
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

subplot('Position',[.775,0.01,.20,.95])
hold off
imagesc(sortrows(ROIclust,j));colormap(Cmap);hold on;axis off
plot([j-.5,j-0.5],[0,size(ROIclust,1)],'k');
plot([j+.5,j+0.5],[0,size(ROIclust,1)],'k');
set(gcf,'Color','w')
% text(j-1,-2,'*','FontSize',20)    

% To better appearance, try organizing data a couple of rounds
foo=ROIclust;
SRV=circshift(SortRowVec,[(j-1),0]);
for k=1:(3*Nlevels)
    foo=sortrows(foo,SRV(k));
end

foo(foo==0)=NaN;
image(sortrows(foo,j));colormap(Cmap);hold on;axis off
plot([j-.5,j-0.5],[0,size(ROIclust,1)],'k');
plot([j+.5,j+0.5],[0,size(ROIclust,1)],'k');
set(gcf,'Color','w')
text(j-0.5,-2,'*','FontSize',20)

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
        fullPath = [fn,sprintf('%0.3f',levels(j)),'.png'];
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