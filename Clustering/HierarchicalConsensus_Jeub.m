function Cons = HierarchicalConsensus_Jeub(S,alpha,N)
% Dependency (see ./ExternalFunction): 
% https://github.com/LJeub/HierarchicalConsensus.git, https://github.com/GenLouvain/GenLouvain.git
% Input:  S assignments nxp where n is nodes and p is partitions
S = double(S);

if ~exist('alpha','var')||isempty(alpha)
    alpha  = 0.05;
end
if ~exist('N','var')||isempty(N)
    N = @(S)localPermModel(S);
end

C = coclassificationMatrix(S);

rng('default');
[Sc,Tree]=hierarchicalConsensus(S,alpha,'CoclassificationMatrix',C,'NullModel',N); % using default settings and alpha = 0.05

[Sall,thresholds] = allPartitions(Sc,Tree);

Cons.allCons = Sall(:,1:end-1)-1; % last one is just single community, we don't need that and because of that we reduce all labels by 1 so the lowest community number start from 1
Cons.allCons=remove_gap(Cons.allCons); % in case there's a gap in the number?
Cons.Sc = Sc;
Cons.Tree = Tree;
Cons.C = C;
Cons.levels = 1:size(Cons.allCons,2);

%% Visualize
Plot_HierachicalConsensus_HSB(Cons,S,C);
end

function newCi = remove_gap(Ci)
all_comms = unique(Ci(:));
newCi= Ci;
for i = 1:length(all_comms)
    newCi(Ci==all_comms(i)) = i;
end
end