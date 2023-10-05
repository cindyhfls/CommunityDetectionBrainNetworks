function [Cons,stats] = Find_Stable_Levels_HSB(stats,partition_centers)
% 2023.09.18 Jiaxin Cindy Tu
%% Remove singleton communities
stats.SortClus=OrgClustMat_HSB(stats.clusters,1); % ouput==0--> junk
subplot(3,1,1);
histogram(stats.SortClus(:),[-1:(1+max(stats.SortClus(:)))])
title(['InfoMap; N=',num2str(length(unique(stats.SortClus(:)))),...
    ' with >=1 members'])
set(gca,'YScale','log')
%% Find pairwise NMI
nclusters = size(stats.SortClus,2);
NMI_dist = zeros(nclusters);
for i = 1:nclusters
    for j = i+1:nclusters
%              tmp = nmi_HSB(double(stats.SortClus(:,i)),double(stats.SortClus(:,j)));
%             NMI_dist(i,j) = 1-tmp.NMI; % somehow this gives NMI >1?
        validid = stats.SortClus(:,i)>0 & stats.SortClus(:,j)>0;
        [~, MIn] = partition_distance(double(stats.SortClus(validid,i)),double(stats.SortClus(validid,j)));
        NMI_dist(i,j) = 1-MIn;
    end
end
NMI_dist = NMI_dist+NMI_dist';

%% reduce the number of partition centers
nmi_th = 0.1; %NMI = 0.9 # this is semi-arbitrary

for k = 1:length(partition_centers)-1
    if k==1
        new_partition_centers = partition_centers(k);
    end
    if NMI_dist( new_partition_centers(end),partition_centers(k+1))>nmi_th
      new_partition_centers =  [new_partition_centers,partition_centers(k+1)];
    end
end

%% Consensus procedure for each

[~,idx] = min(NMI_dist(new_partition_centers,:)); % assign to the closest center

for k = unique(idx)
    Cons.Consensus(:,k) = mode(stats.SortClus(:,idx==k),2);
    Cons.epochs.mean_kden(k)=mean(stats.kdenth(idx==k));
    Cons.epochs.mean_rth(k) = mean(stats.rth(idx==k));
%     Cons.epochs.kden_i(k) = find(stats.kdenth(idx==k),1); %N.B. removed
%     this because the epochs may occasionally not be continous
%      Cons.epochs.kden_f(k) = find(stats.kdenth(idx==k),1,'last');
end

Cons.SortCons = OrgClustMat_HSB(Cons.Consensus,stats.killTH); % remove singleton

%% Visualize NMI in low-D space
xy = cmdscale(NMI_dist,3);

figure('position',[100 100 800 800]);
subplot(2,1,1);hold on;
scatter3(xy(:,1),xy(:,2),xy(:,3),20,stats.kdenth*100,'filled');colormap(jet);c = colorbar('FontSize',15);
c.Label.String = 'edge density';
% h = scatter3(xy(new_partition_centers,1),xy(new_partition_centers,2),xy(new_partition_centers,3),100,'k','p','filled');alpha(0.8);
grid on
title('solution landscape')
set(gca,'YTickLabel',[],'XTickLabel',[],'ZTickLabel',[]);
% legend(h,'solution centers','location','SW');legend('boxoff');
view(45,25);
set(gca,'FontSize',15);
subplot(2,1,2);hold on;
cmap = linspecer(numel(new_partition_centers));
scatter3(xy(:,1),xy(:,2),xy(:,3),20,cmap(idx,:),'filled');
h = scatter3(xy(new_partition_centers,1),xy(new_partition_centers,2),xy(new_partition_centers,3),100,cmap,'p','filled');
grid on
set(gca,'YTickLabel',[],'XTickLabel',[],'ZTickLabel',[]);
legend(h,'solution centers','location','SW');legend('boxoff');
view(45,25);
set(gca,'FontSize',15);
end