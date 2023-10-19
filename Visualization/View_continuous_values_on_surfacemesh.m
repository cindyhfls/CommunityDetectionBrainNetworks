%% Load mesh
function View_continuous_values_on_surfacemesh(vals)

load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr

%% Plot
Lindtrunc= with_without_mw_conversion('Lindtrunc');
Rindtrunc = with_without_mw_conversion('Rindtrunc');
Lindfull= with_without_mw_conversion('Lindfull');
Rindfull = with_without_mw_conversion('Rindfull');
Anat.CtxL.data=zeros(32492,1);Anat.CtxL.data(Lindfull) = vals(Lindtrunc);
Anat.CtxR.data=zeros(32492,1);Anat.CtxR.data(Rindfull-32492) = vals(Rindtrunc);
params.TC = 0;
params.Cmap.P = 'jet';
% cmap = hot(200);
% params.Cmap.P = flipud(cmap(101:200,:));%hot(100);% 'jet';
% cmap = redbluecmap(11);
% params.Cmap.P = flipud(interp1(linspace(0,100,5),cmap(7:11,:),linspace(0,100,100)));
% params.Cmap.P = flipud(interp1(linspace(0,100,10),redbluecmap(10),linspace(0,100,100)));
% cmap = videen(200);
% params.Cmap.P =flipud(cmap(1:100,:));
% params.Cmap.flipP =1;% 'jet';
params.DR=100;
params.Scale=1;
params.Th.P=0;
params.CBar_on =false; % colorbar, the original code had some problem with position
params.PD=1;
params.ctx='inf';           % 'std','inf','vinf'

figure;
subplot(2,1,1);
params.fig_handle = gca;
params.view='lat';       % 'dorsal','post','lat','med'
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
subplot(2,1,2);
params.fig_handle = gca;
params.view='med';       % 'dorsal','post','lat','med'
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
% set(gcf,'color','w');
end