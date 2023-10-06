

%% If loading from saved data
clear;close all;clc;
% filename = './Infomap/Infomap_BCP_220601.mat'
% filename = './Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = './Results/washu120/Gordon/230531/Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup/Results/washu120/Gordon/230531/Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/BCP_Dec_N177/eLABE_Y2_prelim_05062023/230616/Infomap_BCP_Dec_N177_low0.001_step0.001_high0.100_xdist20.mat'
% filename = '/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/WashU120/Gordon/230915/Infomap_WashU120_low0.006_step0.003_high0.150_xdist0.mat';
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
%% Visualize Networks-on-brain and Consensus Edge Density Matrix
Explore_ROI_kden_HSB(foo,CWro.cMap,Anat,params.roi,Cons.epochs.mean_kden);

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
    print([params. outputdir,'/Consensus_Model_',num2str(i)],'-dpng')
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
    temp=Fix_US_Grouping_HSB(Cons.SortConsRO,j); %This code attempts to assign unspecified ROIs to networks, when possible
            temp(string(Nets)=='None'|string(Nets)=='Usp')=0;
    temp=squeeze(Cons.SortConsRO(:,j));
    
    % USp networks with less than 5
%     NnetsI=unique(temp);
%     for nn=1:length(NnetsI)
%         if sum(temp==NnetsI(nn))<5,temp(temp==NnetsI(nn))=0;end
%     end

    temp = Cons.SortConsRO(:,j);
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