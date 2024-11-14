function matOUT = relabel_partitions_HSB(matIN)
% copied from consensus_und, turns the numbers into 1:Nnet
%% reduce to the lowest order
vals = setdiff(unique(matIN),0);
clrs = matIN;
for i=1:size(vals,1)
    clrs(matIN==vals(i))=i;
end
matOUT = clrs;
end