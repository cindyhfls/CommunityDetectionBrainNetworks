function [Cons,stats]=Org_Cons_Org_IMap_Matrix_HSB(stats,cMapTemp)

% This function organizes module labels across edge densities in 
% a matrix of ROI-by-kden output from InfoMap. It then calculates a model
% of the concensus across the matrix. This works best when the steps
% between kden are <=0.01, the smaller the better.
% Then the data is organized again to clean out remaining modules with too
% few members (N-cuttoff is 5 for ROI and 100 for voxels/nodes).

%% Parameters and Initialization
ConsKdenSize=7; % was 5; has been 10, 5, 7 % CT: consecutive edge densities that are treated as a group
f0=figure('Color','w');


%% Organize to >=2
stats.SortClus=OrgClustMat_HSB(stats.clusters,2); % ouput==0--> junk
subplot(3,1,1);
histogram(stats.SortClus(:),[-1:(1+max(stats.SortClus(:)))])
title(['InfoMap; N=',num2str(length(unique(stats.SortClus(:)))),...
    ' with >=2 members'])
set(gca,'YScale','log')

%% Consensus
if ~exist('cMapTemp','var')
    cMapTemp=jet(max(stats.SortClus(:)));
end
[Cons]=Consensus_Auto_HSB(stats,cMapTemp,ConsKdenSize);
figure(f0)
subplot(3,1,2);
histogram(Cons.Consensus(:),[-1:(1+max(Cons.Consensus(:)))])
title(['Consensus; N=',num2str(length(unique(Cons.Consensus(:)))),...
    ' unique modules'])
set(gca,'YScale','log')

%% Organize to >=5
Cons.SortCons=OrgClustMat_HSB(Cons.Consensus,stats.params.killTH);
subplot(3,1,3);
histogram(Cons.SortCons(:),[-1:(1+max(Cons.SortCons(:)))])
title(['Consensus; N=',num2str(length(unique(Cons.SortCons(:)))),...
    ' unique modules with >=',num2str(stats.params.killTH),' members'])
set(gca,'YScale','log')

for j=1:size(Cons.SortCons,2)
Cons.epochs.mean_kden(j)=mean(stats.kdenth(...
    Cons.epochs.kden_i(j):Cons.epochs.kden_f(j)));
Cons.epochs.mean_rth(j)=mean(stats.rth(...
    Cons.epochs.kden_i(j):Cons.epochs.kden_f(j)));
end

%% Stats of consensus models
Cons.stats=Matrix_metrics_HSB(Cons.SortCons,FisherZ2R_HSB(stats.MuMat),...
    Cons.epochs.mean_rth,stats.params.binary);

Ncons=size(Cons.SortCons,2);

%% Plot
figure('position',[100 100 800 800]);
subplot(2,2,1);
plot([1:Ncons],Cons.stats.modularity,'bo-','LineWidth',2);
axis([0,Ncons+1,0,1.1]);grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Modularity'},'Location','southeast');

subplot(2,2,2);
plot([1:Ncons],Cons.stats.NnBc,'ro-','LineWidth',2);
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Proportion of nodes in the largest component'},'Location','southeast');

subplot(2,2,3);
plot([1:Ncons],Cons.stats.AvgSil,'ko-','LineWidth',2);
grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Average Silhouette'},'Location','southeast');

subplot(2,2,4);
plot([1:Ncons],Cons.stats.num_communities,'ko-','LineWidth',2);
grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Number of Communities'},'Location','southeast');


