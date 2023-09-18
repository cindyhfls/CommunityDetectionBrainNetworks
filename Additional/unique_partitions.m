function [ciu,unique_counts] = unique_partitions(ci)
% copied and modified from consensus_und
ci = relabel_partitions(ci);
ciu = [];
count = 0;
c = 1:size(ci,2);
unique_counts = [];
while ~isempty(ci)
    count = count + 1;
    tgt = ci(:,1);
    ciu = [ciu,tgt];                %#ok<AGROW>
    dff = sum(abs(bsxfun(@minus,ci,tgt))) == 0;
    ci(:,dff) = [];
    c(dff) = [];
    unique_counts(count) = sum(dff);
end