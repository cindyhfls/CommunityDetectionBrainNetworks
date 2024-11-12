function matIN = remove_singleton(matIN,killTH)
% this scripts removes singleton communities or communities with size
% smaller than killTH
if ~exist('killTH','var')||isempty(killTH)
    killTH = 2;
end
%% Set grouping value == 0 if less than threshold membership
[~,Nkden]=size(matIN);
spread=zeros(Nkden,1);  % start with end col with fewest clusters
for j=1:Nkden
    vals=unique(squeeze(matIN(:,j)));
    spread(j)=length(vals);
    for k=1:spread(j)
        idx=(matIN(:,j)==vals(k));
        if sum(idx)<killTH,matIN(idx,j)=0;end
    end
end
%% reduce to the lowest order
vals = setdiff(unique(matIN),0);
clrs = matIN;
for i=1:size(vals,1)
    clrs(matIN==vals(i))=i;
end
matIN = clrs;
end