function f =View_prop_Colors_transparent(thisprop,Parcels)

load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr
%%
cmap = flipud(bone(100));
f = figure('position',[100 100 385 275]);
ax1 = subplot(2,1,1);
set(ax1,'Position',[0.01 0.5,0.91,0.5]);
plot_parcels_by_values(thisprop,Anat,'lat',Parcels,[0,1],cmap);
ax2 = subplot(2,1,2);
set(ax2,'Position',[0.01,0.05,0.91,0.5]);
plot_parcels_by_values(thisprop,Anat,'med',Parcels,[0,1],cmap);
colormap(cmap);caxis([0,1]);
colorbar('position',[0.93, 0.1, 0.02, 0.85],'Ticks',[0,1]);
set(gca,'FontSize',18);
end