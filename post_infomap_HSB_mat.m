

%% If loading from saved data
clear;close all;clc;
% filename = './Infomap/Infomap_BCP_220601.mat'
% filename = './Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = './Results/washu120/Gordon/230531/Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup/Results/washu120/Gordon/230531/Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/BCP_Dec_N177/eLABE_Y2_prelim_05062023/230616/Infomap_BCP_Dec_N177_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/WashU120/Gordon/230904/Infomap_WashU120_low0.006_step0.003_high0.150_xdist20.mat';
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/Gordon/2310904/Infomap_eLABE_Y2_N113_low0.010_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/eLABE_Y2_prelim_072023_0.75/230927/Infomap_eLABE_Y2_N113_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/eLABE_Y2_prelim_072023_0.75/230904/Infomap_eLABE_Y2_N113_low0.006_step0.001_high0.150_xdist0.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/WashU120/Gordon/231129/Infomap_WashU120_low0.010_step0.001_high0.200_xdist20.mat';
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/Gordon/231011/Infomap_eLABE_Y2_N113_low0.010_step0.001_high0.300_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/WashU120/Gordon/231016/Infomap_WashU120_low0.006_step0.001_high0.200_xdist20.mat';
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N92_healthyterm/Tu_342/231106/Infomap_eLABE_Y2_N92_healthyterm_low0.010_step0.001_high0.200_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/Tu_342/231016/Infomap_eLABE_Y2_N113_low0.010_step0.001_high0.200_xdist20.mat'
filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N92_healthyterm/Tu_326/240522/Infomap_eLABE_Y2_N92_healthyterm_low0.003_step0.003_high0.200_xdist30.mat'


load(filename)
% stats = stats{1}; % at sometime in 202205 I changed the Infomap stats results to a cell format so I can save multiple results
params = stats.params

Nroi = size(params.roi,1)

% if loading mat again after clearing the space comment the two lines below out
if ~isfield(stats,'MuMat')||isempty(stats.MuMat)
    tmp = smartload(stats.params.zmatfile);
    stats.MuMat = mean(tmp,3);
end
% figdir = fullfile('./Figures',params.IMap_fn);

parcel_name =params.parcel_name%'eLABE_Y2_prelim_072023_0.75'%'Gordon'% params.parcel_name
load(['Parcels_',parcel_name,'.mat'],'Parcels');
%% Sorting the solutions sequentially % N.B. original code by J. Powers occassionally changes the community assignment so use the one from genlouvain
minsize = 5;
stats.SortClus = postprocess_ordinal_multilayer(stats.clusters);
stats.SortClus = remove_singleton(stats.SortClus,minsize);
stats.SortClus = remove_few_column_clusters(stats.SortClus);
%% Sort all densities and assign colors

nameoption = 3% 1: automatic, 3: using template
% templatepath = '/data/wheelock/data1/parcellations/Myers-Labonte2024/elabe_n131_parcels_height0.5.dlabel.nii'
templatepath  ='Myers-Labonte2024_23Networks.dlabel.nii'%'Gordon2017_17Networks.dlabel.nii'% 'Tu_eLABE_Y2_22Networks.nii'
parcelpath ='/data/wheelock/data1/parcellations/InfantParcellation_Tu/May2024/eLABE_Y2_N92_healthyterm_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.65_heightperc_0.9_minsize_15_relabelled.dlabel.nii'
% parcelpath ='/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/eLABE_Y2_N113_atleast600frames_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.75_nogap_minsize_15_relabelled.dlabel.nii';
% parcelpath ='/data/wheelock/data1/parcellations/InfantParcellation_Tu/Oct2023/eLABE_Y2_N113_atleast600frames_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.65_heightperc_0.9_minsize_15_relabelled_N342.dlabel.nii';
% parcelpath = '/data/wheelock/data1/parcellations/333parcels/Parcels_LR.dtseries.nii'
% [CWro,stats.SortClusRO] = assign_network_colors(stats.SortClus,nameoption); % currently using Gordon 13 network colors as default
[CWro,stats.SortClusRO] = assign_network_colors(stats.SortClus,nameoption,templatepath,parcelpath);

%% Plot single column
for thresh = [0.25,1.25,1.75,2.75,4,5.75,10,16,19.25]
    icol = find(round(stats.kdenth*100,2)==thresh)
    plot_network_assignment_parcel_key(Parcels, stats.SortClusRO(:,icol),CWro.cMap,CWro.Nets)
    print(fullfile(params.outputdir,sprintf('kden_%1.4f.tif',stats.kdenth(icol))),'-dtiff','-r300');
    close all
end
%% (optional) Viewing and Manual edit of specific networks/Save out those networks
for iNet =1:length(CWro.Nets)
    Edit_NetworkColors(stats.SortClusRO,CWro,iNet,Parcels);
    minthresh = min(stats.kdenth(any(stats.SortClusRO==iNet)))*100;
    maxthresh = max(stats.kdenth(any(stats.SortClusRO==iNet)))*100;
    text(0.2,1.2,sprintf('%1.2f%%-%1.2f%%',minthresh,maxthresh),'Units','Normalized');
    print(fullfile(params.outputdir,sprintf('Networks%02d',iNet)),'-dtiff','-r300');
    %     pause;
    close all;
end
% customcolor = [];% fill in if you want to change it
% customcolor = distinguishable_colors(1,CWro.cMap); % set color to
% one
% CWro=Edit_NetworkColors(newstats.SortSortClus,CWro,iNet,Parcels,customcolor);
%% Visualize Networks-on-brain and Consensus Edge Density Matrix

% Explore_ROI_levels_HSB(foo,CWro.cMap,Anat,params.roi,Cons.epochs.mean_kden);
init_level = find(abs(stats.kdenth-0.0275)<1e-4,1); % start at a specific level
% init_level = 1 % start from first column
Explore_parcel_levels_HSB(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'kden'),init_level);

%% Plot a legend for the Networks
N = length(CWro.Nets)
% figure('Units','inches','position',[10 10 2,3]);%[10 10 5 2]
figure('Units','inches','position',[10 10 5,3]);%[10 10 5 2]
h = gscatter(ones(1,N),ones(1,N),CWro.Nets,CWro.cMap,'s',50);
for i = 1:N
    set(h(i),'Color','k','MarkerFaceColor',CWro.cMap(i,:));
end
legend(CWro.Nets,'interpreter','none','FontSize',10,'location','best','Orientation','horizontal','NumColumns',2);
legend('boxoff')
xlim([10,11]);
axis('off')
print(fullfile(params.outputdir,'networks_Legend'),'-dtiff','-r300');
%% Make video
Make_parcel_kden_Video(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,strrep(params.IMap_fn,'.mat','')))

%% Consensus Simple (adding higher threshold assignments to lower thresholds)
mincol = 1; minsize = 5;
consensusmap = Consensus_infomap_simple_HSB(stats.SortClusRO,mincol,minsize);
Explore_parcel_levels_HSB(consensusmap,CWro.cMap,Parcels,0,fullfile(params.outputdir,'consensus_simple'));

%% Save the consensus as text and dlabel file for manual editing
parcelcifti = cifti_read(parcelpath);
Nparcels = size(stats.clusters,1);
for k = 1:Nparcels
    currcolor = CWro.cMap(consensusmap(k),:);
    parcelcifti.diminfo{1,2}.maps.table(k+1).name = sprintf('%i: %s',k,CWro.Nets{consensusmap(k)});
    parcelcifti.diminfo{1,2}.maps.table(k+1).rgba = [currcolor';1];
end
cifti_write(parcelcifti,fullfile(params.outputdir,'consensus_simple.dlabel.nii'));
writematrix(consensusmap,fullfile(params.outputdir,'consensus_simple.txt'))

%% Now find the consensus
[Cons] = Find_Stable_Levels_HSB(stats); % consensus by finding
% consecutive stable levels from the NMI, Eggebrecht et al. 2017 Cerebral
% Cortex but determin the levels in a data driven way rather than
% pre-defined
Cons = Cons_metrics_HSB(Cons,stats); % get some stats for the consensus and plot the figure

save(fullfile(params.outputdir,params.IMap_fn),'Cons','-append')

Explore_parcel_levels_HSB(Cons.SortCons,CWro.cMap,Parcels,Cons.kdenth,fullfile(stats.params.outputdir,'consensus'));

%%

S = horzcat(stats.allpartitions{:});

[Cons] = HierarchicalConsensus_Jeub(S,0.05,@(S)permModel(S)); % consensus with Jeub et al. 2018 Scientific Reports
Plot_HierachicalConsensus_HSB(Cons,S,Cons.C);
print(gcf,fullfile(params.outputdir,strrep(params.IMap_fn,'.mat','_CoassignmentMatrix')),'-dtiff','-r300');

tmp = postprocess_ordinal_multilayer(fliplr(Cons.SortCons));
Cons.SortCons = fliplr(tmp);
Cons.SortCons= remove_singleton(Cons.SortCons,2);

save(fullfile(params.outputdir,params.IMap_fn),'Cons','-append')

% Assign color to consensus
nameoption = 3;% 1: automatic, 3: using template
templatepath  ='Gordon2017_17Networks.dlabel.nii';% 'Tu_eLABE_Y2_22Networks.nii'
% parcelpath ='/data/wheelock/data1/parcellations/InfantParcellation_Tu/Oct2023/eLABE_Y2_N113_atleast600frames_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.65_heightperc_0.9_minsize_15_relabelled_N342.dlabel.nii';
parcelpath = '/data/wheelock/data1/parcellations/333parcels/Parcels_LR.dtseries.nii'
[CWro,Cons.SortConsRO] = assign_network_colors(Cons.SortCons,3,templatepath,parcelpath);

%%
warning('off');
Explore_parcel_levels_HSB(Cons.SortConsRO,CWro.cMap,Parcels,Cons.levels,fullfile(stats.params.outputdir,'consensus'));
%  Explore_parcel_levels_HSB(Cons.modeCons,CWro.cMap,Parcels,Cons.mean_kdenth);
close all;

Cons = Cons_metrics_HSB(Cons,stats); % get some stats for the consensus and plot the figure
print(gcf,fullfile(params.outputdir,strrep(params.IMap_fn,'.mat','_Consensus_metrics')),'-dtiff','-r300');

%%
Zr = NaN(stats.params.repeats);
for i = 1:stats.params.repeats
    for j = i+1:stats.params.repeats
        Zr(j,i) = zrand(stats.allpartitions{ithre}(:,i),stats.allpartitions{ithre}(:,j));
    end
end

%% Plot Fc matrix
icons = 4
key = Cons.SortConsRO(:,icons);key(key==0) = length(CWro.Nets)+1;
CWro.cMap(length(CWro.Nets)+1,:)=[0.5,0.5,0.5];
% first sort in the order that goes in
[~,sortid] = sort(key);

figure;
% imagesc(stats.MuMat(sortid,sortid));
Matrix_Org3(stats.MuMat(sortid,sortid),repmat(key(sortid),1,2),10,[-1,1],CWro.cMap,0,jet)
c = colorbar;
c.Label.String = 'z(r)'
set(gca,'FontSize',15);
print(fullfile(params.outputdir,['ConsensusMatrixPlot_',num2str(icons)]),'-dtiff','-r300');
%% Plot spring-embedded plot?
stats.MuMat;
G = graph(thresholded_matrix,'upper');% sometimes the matrix is not symmetric? precision problem?
Lwidths = 1*G.Edges.Weight/max(G.Edges.Weight)*0.05;
figure;
h = plot(G,'ko-','layout','force','UseGravity',true,'NodeCData',stats.SortClusRO(:,1),'NodeColor','flat','MarkerSize',2,'LineWidth',Lwidths);
colormap(CWro.cMap)

% or use this
% spring_embedding_func_easy_crossthresh(corrmat,assignments,Kc,L0,distances,xdist,thresholds,outname)
% Kc: try 1
% L0:try 25
return
%% Simple consensus
minsize =5;
lowestcol = 1%find(abs(stats.kdenth-single(0.01))<10E-5);
highestscol =21%find(abs(stats.kdenth-single(0.10))<10E-5)%size(stats.clusters,2);
stats.SortClus =OrgClustMat_HSB(stats.clusters,minsize,0); % last argument = 1 for reverse ordering
consensusmap = Consensus_infomap_simple_CT_mod(stats.SortClus,lowestcol,minsize);

Cons.SortClus = OrgClustMat_HSB(consensusmap ,minsize,0); 
% templatepath = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup/Templates/Tu_eLABE_Y2_22Networks.nii';%Laumann2015_12Networks.dlabel.nii'
% parcelpath = '/data/wheelock/data1/parcellations/333parcels/Parcels_LR.dtseries.nii'
% parcelpath ='/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/eLABE_Y2_N113_atleast600frames_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.75_nogap_minsize_15_relabelled.dlabel.nii';
[CWro,Cons] = assign_network_colors(Cons,3)%,templatepath,parcelpath) % currently using Gordon 13 network colors as default

% [CWro,Cons] = assign_network_colors(Cons,3) % currently using Gordon 13 network colors as default
pause(0.1);close all

parcel_name =params.parcel_name%'eLABE_Y2_prelim_072023_0.75'%'Gordon'% params.parcel_name
load(['Parcels_',parcel_name,'.mat'],'Parcels');
figure('position',[100 100 400 300]);
for i = 1:size(Cons.SortClusRO,2) 
    key =Cons.SortClusRO(:,i);
    plot_network_assignment_parcel_key(Parcels, key,CWro.cMap,CWro.Nets,0)  
    text(0.72,0,'simple consensus','Units','Normalized')
%     print([params.outputdir,'/Consensus_Model_SimpleConsensus'],'-dtiff','-r300')
%     pause;
%     clf
end

cMap=CWro.cMap;
Nets=CWro.Nets;

temp = Cons.SortClusRO;
if any(temp==0)
    if sum(contains(Nets,{'None','Usp'}))
        noneidx = find((string(Nets)=="None")|(string(Nets)=="USp"));
        temp(temp==0) = noneidx;
    else
    temp(temp==0)=size(cMap,1)+1; % Unspecified network became the last network
    cMap=cat(1,cMap,[0.25,0.25,0.25]);% gray for USp
    Nets=cat(1,Nets,'None');
    end
end
keep=unique(temp)';

% Put together IM structure
clear IM
IM.name=['IM_',params.IMap_fn,'_Consesus_model_simple'];
IM.cMap=cMap(keep,:);
IM.Nets=Nets(keep);
IM.ROIxyz=params.roi;
IM.key=[[1:Nroi]',zeros(Nroi,1)];
[IM.key(:,2),IM.order]=sort(IM_Remove_Naming_Gaps_HSB(temp));
IM.ROIxyz=IM.ROIxyz(IM.order,:);
IM=Org_IM_DVLR_HSB(IM);

% Visualize
figure;
noneidx = find((string(IM.Nets)=="None")|(string(IM.Nets)=="USp"));
if ~isempty(noneidx)
    keepnets =IM.key(:,2)~=noneidx;
else
    keepnets = true(size(IM.key(:,2)));
end
M = ones(max(IM.key(:,2)));M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
zmat = stats.MuMat(IM.order,IM.order);
Matrix_Org3(zmat(keepnets,keepnets),...
    IM.key(keepnets,:),10,[-0.6,0.6],IM.cMap,0);

% Matrix_Org3_HSB(stats.MuMat(IM.order,IM.order),IM.key,10,[-0.6,0.6],IM.cMap,0);
% title(sprintf('%s, %i Networks',[strrep(IM.name,'_',' ')],length(IM.Nets)),'Color','w');

D = calc_correlationdist(stats.MuMat(IM.order,IM.order));
s = silhouette_mod(IM.key(keepnets,2),D(keepnets,keepnets),M);
text(0.3,-0.05,sprintf('avg SI = %2.3f',mean(s)),'Units','normalized','FontWeight','Bold');
% print(fullfile(params.outputdir,'Heatmap_SimpleConsensus'),'-dtiff','-r300')
%% plot the versatility
if params.repeats_consensus
    figure('position',[100 100 800 400]);
    yyaxis right;
    errorbar(stats.kdenth,stats.avgVersatility,stats.stdVersatility/sqrt(size(stats.clusters,1)));
    ylabel('Versatility');
    xlabel('kdenth');
    legend('SEM')
    yyaxis left;
    plot(stats.kdenth,stats.metrics.non_singleton);
    ylabel('# communities');
end
%% Plot other metrics
figure('position',[100 100 800 400]);hold on;
yyaxis right;
errorbar(stats.kdenth,stats.metrics.AvgSil,stats.metrics.StdSil/sqrt(size(stats.clusters,1)))
ylabel('Silhouette Index');
yyaxis left;
plot(stats.kdenth,stats.metrics.non_singleton);
ylabel('# communities');
%% Plot DB and CH
figure;
subplot(2,1,1);
plot(stats.kdenth,stats.metrics.DB);
title('Davies-Bouldin Criterion');
subplot(2,1,2);
plot(stats.kdenth,stats.metrics.CH);
title('CalinskiHarabasz Criterion')

return;
%% Identify the local minimum of versatility as a typical solution
% find(islocalmin(stats.avgVersatility,'FlatSelection','first'))
% find(islocalmin(stats.metrics.DB,'FlatSelection','first'))
% find(islocalmax(stats.metrics.CH,'FlatSelection','first'))
% find(islocalmax(stats.metrics.AvgSil,'FlatSelection','first'))
partition_centers =find(islocalmin(stats.avgVersatility,'FlatSelection','first')|...
    islocalmax(stats.metrics.AvgSil,'FlatSelection','first')|...
    islocalmin(stats.metrics.DB,'FlatSelection','first')|...
    islocalmax(stats.metrics.CH,'FlatSelection','first'))
%% Run Consensus Procedure to reduce the number of networks -> planning to replace this
% uses normalized mutual information across thresholds specified in params

% Adam's original method (finding stable neighboring threshold with at least x consecutive)
[Cons,stats]=Org_Cons_Org_IMap_Matrix_HSB(stats,[],3); % this is Adam's original workflow to find stable groups with high NMI in the neighboring thresholds, I will update that with a different function
% Updated to calculate pairwise NMI
[Cons,stats] = Find_Stable_Levels_HSB(stats,partition_centers); % consensus by finding stable levels from the 

Cons = Cons_stats_HSB(Cons,stats); % get some stats for the consensus and plot the figure
% stats.SortedStats = Matrix_metrics_HSB(stats.SortClus,stats.MuMat,stats.rth,stats.params.binary,stats.params.type,stats.kdenth);

% save(filename,'Cons','-append'); % save infomap output to matrix

%% load MNI mesh
load('MNI_coord_meshes_32k.mat','MNIl','MNIr');
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
clear MNIl MNIr

%%
nameoption = 3
[CWro,Cons] = assign_network_colors(Cons,nameoption) % currently using Gordon 13 network colors as default
foo=Cons.SortConsRO;

%% (Alternatively) Visualize the parcels
parcel_name ='Gordon'% params.parcel_name
load(['Parcels_',parcel_name,'.mat'],'Parcels');
figure('position',[100 100 400 300]);

for i = 1:size(foo,2)
    key = foo(:,i);
    plot_network_assignment_parcel_key(Parcels, key,CWro.cMap,CWro.Nets)  
  text(0.72,0,sprintf('Avg density = %2.2f%%',Cons.epochs.mean_kden(i)*100),'Units','Normalized')
%     plot_network_assignment_parcel_key(Parcels, key,CWro.cMap,CWro.Nets)
%     text(0.1,1.5,sprintf('Avg density = %2.2f%%',Cons.epochs.mean_kden(i)*100),'Units','Normalized')
    print([params. outputdir,'/Consensus_Model_',num2str(i)],'-dtiff')
%     pause;
    clf
end
%% Generate Infomap (IM) Structure for Viable Edge Density Ranges
% Viable = edge densities in which connectivity >80% (see figure output
% from Org_Cons_Org_Imap_Matrix)
% IM structures are used during Enrichment to organize ROI into networks
% This set of codes visualizes all possible IM options to choose from

% load FC matrix
if strcmp(stats.params.format,'mat')
    stats.MuMat = smartload(stats.params.zmatfile); %(parcel-wise data in .mat)
elseif strcmp(stats.params.format,'cifti')
    tmp = cifti_read(stats.params.zmatfile);stats.MuMat = tmp.cdata; %(vertex-wise datain cifti format)
end

toIM=[1:size(Cons.SortCons,2)];
for j=toIM % Auto out of IM for each Cons model
    
    %     if (Cons.stats.NnBc(j)>0.9) && (Cons.stats.kave(j)>log(Nroi))
    % remove number of nodes in largest component <= 90%? and mean degree
    % <log(Nroi)? degree is calculated with number of non-zero weight connections
    
    
    % Turn into function? inputs: Cons, CWro, IM name, stats
    
    cMap=CWro.cMap;
    Nets=CWro.Nets;
    
    % Add a way to fix USp?
    temp=Fix_US_Grouping_HSB(Cons.SortClusRO,j); %This code attempts to assign unspecified ROIs to networks, when possible
            temp(string(Nets)=='None'|string(Nets)=='Usp')=0;
    temp=squeeze(Cons.SortClusRO(:,j));
    
    % USp networks with less than 5
%     NnetsI=unique(temp);
%     for nn=1:length(NnetsI)
%         if sum(temp==NnetsI(nn))<5,temp(temp==NnetsI(nn))=0;end
%     end

    temp = Cons.SortClusRO(:,j);
    if any(temp==0)
        temp(temp==0)=size(cMap,1)+1; % Unspecified network became the last network
        cMap=cat(1,cMap,[0.25,0.25,0.25]);% gray for USp
        %             cMap = cat(1,cMap,[1,1,0.8]); % a very light yellow for USp
        Nets=cat(1,Nets,'None');
    end
    keep=unique(temp)';
    
    % Put together IM structure
    clear IM
    IM.name=['IM_',params.IMap_fn,'_Consesus_model_',num2str(j)];
    IM.cMap=cMap(keep,:);
    IM.Nets=Nets(keep);
    IM.ROIxyz=params.roi;
    IM.key=[[1:Nroi]',zeros(Nroi,1)];
    [IM.key(:,2),IM.order]=sort(IM_Remove_Naming_Gaps_HSB(temp));
    IM.ROIxyz=IM.ROIxyz(IM.order,:);
    IM=Org_IM_DVLR_HSB(IM);    
    
    % Visualize
    figure('Color','k','Units','Normalized','Position',[0.35,0.30,0.35,0.61]);
    subplot(4,1,[1:3])
    Matrix_Org3_HSB(stats.MuMat(IM.order,IM.order),IM.key,10,[-0.6,0.6],IM.cMap,0);
    title(sprintf('%s,kden=%0.2f%% %i Networks',[strrep(IM.name,'_',' ')],Cons.epochs.mean_kden(j)*100,length(IM.Nets)),'Color','w');
    
%     title([{[strrep(IM.name,'_',' ')]};{['kden=',...
%         num2str(Cons.epochs.mean_kden(j)*100,'%0.2f'),'%, ',...
%         num2str(stats.kdenth(Cons.epochs.kden_i(j))*100,'%0.2f'),...
%         '-',num2str(stats.kdenth(Cons.epochs.kden_f(j))*100,'%0.2f'),'%, ',...
%         num2str(length(IM.Nets)),' Networks']}],'Color','w')
    c=colorbar;c.Ticks=[-0.6,0,0.6];c.TickLabels={'-0.6','0','0.6'};
    set(c,'Color','w')
    subplot(4,1,4)
    histogram(IM.key(:,2));title(strrep(IM.name,'_',' '))
    set(gca,'XTick',[1:max(IM.key(:,2))],'XTickLabel',...
        IM.Nets,'Color','k',...
        'XColor','w','YColor','w');
    xtickangle(45);
    ylabel('Nrois','Color','w')
    xlim([0,max(IM.key(:,2)+1)]);
    set(gcf,'InvertHardCopy','off');
%     print(gcf,fullfile(params.outputdir,sprintf('%s_heatmap',IM.name)),'-dtiff')
    
    % if function, end here
    
    % Save the IM file
%     save(fullfile(params.outputdir,[IM.name,'.mat']),'IM');
    
    %     end
end

%% Visualize IM Model on Brain with Network Names and Colors
params.radius = 4;
Anat.alpha = 1;
% IM.cMap(end,:) = [0.5,0.5,0.5];% set unassigned to gray

figure; % this shows the sorted FC
Matrix_Org3_HSB(stats.MuMat(IM.order,IM.order),...
    IM.key,10,[-0.3,0.3],IM.cMap,1); % mean

figure; % this shows the center of the parcels in a sphere with radii params.radius
Anat.ctx='std';View_ROI_Modules(IM,Anat,IM.ROIxyz,params);
%% Visualize some stats
figure;
subplot(2,2,1);
plot(stats.kdenth*100,stats.metrics.non_singleton);
subplot(2,2,2);
plot(stats.kdenth*100,stats.metrics.Nc);
subplot(2,2,3);
plot(stats.kdenth*100,stats.metrics.Cdns);
subplot(2,2,4);
plot(stats.kdenth*100,stats.metrics.AvgSil);