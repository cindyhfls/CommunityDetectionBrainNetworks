function prop = calc_network_prop(Clust,G1)
%% Calculate the proportion of subject with network j at each parcel
% G1=setdiff(unique(Clust(:)),0); % because sometimes they may not have all
% the clusters so we manually input it
prop = NaN(size(Clust,1),length(G1));
for j = 1:length(G1)
    prop(:,j) = mean(Clust==G1(j),2);
end
prop(prop==0) = NaN; % set to remove
end