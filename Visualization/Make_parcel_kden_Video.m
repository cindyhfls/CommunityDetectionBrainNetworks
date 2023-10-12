function Make_parcel_kden_Video(ROIclust,Cmap,Parcels,kden,fn)
%
% This function generates a movie of ROI Sortings for
% various kden, scanning through the kden chosen.
if any(ROIclust(:)==0)
    Cmap=cat(1,Cmap,[0.5,0.5,0.5]);
    ROIclust(ROIclust==0)=size(Cmap,1);
end
%% Set up paramters
[Nroi,Nkden]=size(ROIclust);
params.Cmap.P=Cmap;%IM.cMap;jet(nNet)
params.TC=1;
params.ctx='inf';         % also, 'std','inf','vinf'

writerObj = VideoWriter(fn,'Motion JPEG AVI');%VideoWriter(fn,'MPEG-4');
writerObj.FrameRate = 4;
writerObj.Quality = 100;

%% Display figures and write movie
open(writerObj);
f=figure('Position',[310,530,1550,570]);
for j=1:Nkden
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
%% Write to video
% pause(2)
drawnow;
frame = getframe(f);
writeVideo(writerObj,frame);
end
close(writerObj);
