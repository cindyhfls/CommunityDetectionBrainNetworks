function Plot_HierachicalConsensus_HSB(Cons,S,C)
% input: Cons - structure containing Sc (finest level of partition) and
% Tree (hierarchical tree indicating how to merge clusters in Sc to
% reconstruct coarser clusters, see hierchicalConsensus.m
%
Sc = Cons.Sc;
Tree = Cons.Tree;
%% Visualization
if ~exist('C','Var')||isempty(C)
    C = coclassificationMatrix(S);
end

figure;
consensusPlot(C,Sc,Tree);

end