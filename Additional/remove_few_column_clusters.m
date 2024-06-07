function matIN = remove_few_column_clusters(matIN,killTH)
% this scripts removes communities that only exist in one column (e.g. edge
% threshold), this code is a separate script from remove_singleton just in
% case there are reasons to keep the clusters that exist in a few threshold
% (e.g. only have one column)

if ~exist('killTH','var')||isempty(killTH)
    killTH = 2;
end
%% Count the number of thresholds
allnets=reshape(setdiff(unique(matIN(:)),0),1,[]);
for inet = allnets
    if length(find(any(matIN==inet)))<killTH
        matIN(matIN==inet)=0;
    end
end
%% reduce to the lowest order
vals = setdiff(unique(matIN),0);
clrs = matIN;
for i=1:size(vals,1)
    clrs(matIN==vals(i))=i;
end
% clrs=single(clrs); % Changing data type to single seems to cause more trouble downstream than the
% slight saving in storage
matIN = clrs;
end