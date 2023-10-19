function [Cons,stats] = Find_Stable_Levels_HSB(stats,CWro,Parcels)
% 2023.09.18 Jiaxin Cindy Tu
%% Find pairwise NMI
nclusters = size(stats.clusters,2);
nroi = size(stats.clusters,1);
NMI_dist = zeros(nclusters);
for i = 1:nclusters
    for j = i+1:nclusters
        [~,MIn] = partition_distance(double(stats.clusters(:,i)),double(stats.clusters(:,j))); % use the origional cluster assignment not the cleaned one
        NMI_dist(i,j) = 1-MIn;
    end
end
NMI_dist = NMI_dist+NMI_dist';
%% Plot NMI distance
d_before = arrayfun(@(i) NMI_dist(i,i-1),2:length(NMI_dist)-1);
borders = find(isoutlier(d_before,'median')); % outliers identified as three median absolute deviations away from the median

borders = unique([1,borders,length(NMI_dist)+1]);
epochs = cell(length(borders)-1,1);
keepidx =true(size(epochs));
for i = 1:length(borders)-1
    epochs{i} = borders(i):borders(i+1)-1;
    if length(epochs{i})==1 % only one threshold
        keepidx(i) = false;
    end
end

% remove the epochs with only one threshold
epochs= epochs(keepidx);
borders = borders(keepidx);

grpclrs = distinguishable_colors(length(epochs));
figure('position',[100 100 1200 400]);
subplot(1,3,1:2);hold on;
xx = 1:length(NMI_dist);
plot(xx(2:length(NMI_dist)-1),1-d_before,'k-')
yl = ylim;
for ii = 1:size(grpclrs,1)
    patch([epochs{ii}(1),epochs{ii}(end),epochs{ii}(end),epochs{ii}(1)],[yl(1),yl(1),yl(end),yl(end)],grpclrs(ii,:),'LineStyle','None');alpha(0.4);
end
% arrayfun(@(ii)vline(ii,'LineStyle','-','Color','r'),xx(borders));
ylabel('NMI');
xlabel('threshold number');
title('NMI to the solution in an earlier threshold');

subplot(1,3,3)
imagesc(NMI_dist);
arrayfun(@(ii)vline(ii),borders);
arrayfun(@(ii)hline(ii),borders);
axis off
axis square
colorbar;
title('NMI distance');
print(fullfile(stats.params.outputdir,'NMIgroups.png'),'-dpng');
%% Examine each epoch

% warning('off');
% for i = 1:length(borders)-1
%     Explore_parcel_kden_HSB(stats.SortClusRO(:,epochs{i}),CWro.cMap,Parcels,stats.kdenth(epochs{i}));
%     close all;
% end

%% Consensus by similarity to all clusters in the epoch
Cons.epochs = epochs;
Cons.SortCons =NaN(nroi,length(epochs));
for i = 1:length(epochs)
    [~,cons_id] = min(sum(NMI_dist(epochs{i},epochs{i}))/(length(epochs{i})-1));
    Cons.centers(i,1) = single(epochs{i}(cons_id));
    Cons.mean_kdenth(i,1) = stats.kdenth(Cons.centers(i,1));
    Cons.mean_rth(i,1) = stats.rth(Cons.centers(i,1));
    Cons.SortCons(:,i) =single(stats.SortClusRO(:,Cons.centers(i,1))); % taking the most representative one
%     Cons.modeCons(:,i) = mode(stats.SortClusRO(:,epochs{i}),2);
end

warning('off');
Explore_parcel_kden_HSB(Cons.SortCons,CWro.cMap,Parcels,Cons.mean_kdenth,fullfile(stats.params.outputdir,'consensus'));
%  Explore_parcel_kden_HSB(Cons.modeCons,CWro.cMap,Parcels,Cons.mean_kdenth);
close all;


%% Visualize NMI in low-D space
xy = cmdscale(NMI_dist,3);
new_partition_centers = Cons.centers;
figure('position',[100 100 800 800]);
subplot(2,1,1);hold on;
scatter3(xy(:,1),xy(:,2),xy(:,3),20,stats.kdenth*100,'filled');colormap(jet);c = colorbar('FontSize',15);
c.Label.String = 'edge density (%)';
% h = scatter3(xy(new_partition_centers,1),xy(new_partition_centers,2),xy(new_partition_centers,3),100,'k','p','filled');alpha(0.8);
grid on
title('solution landscape')
set(gca,'YTickLabel',[],'XTickLabel',[],'ZTickLabel',[]);
% legend(h,'solution centers','location','SW');legend('boxoff');
view(45,25);
set(gca,'FontSize',15);

subplot(2,1,2);hold on;
for i = 1:length(epochs)
    idx = epochs{i};
    scatter3(xy(idx,1),xy(idx,2),xy(idx,3),20,grpclrs(i,:),'filled');
end
h = scatter3(xy(new_partition_centers,1),xy(new_partition_centers,2),xy(new_partition_centers,3),100,'k','p','filled');

grid on
set(gca,'YTickLabel',[],'XTickLabel',[],'ZTickLabel',[]);
legend(h,'solution centers','location','SW');legend('boxoff');
view(45,25);
set(gca,'FontSize',15);
print(fullfile(stats.params.outputdir,'solution_landscape.png'),'-dpng');
end