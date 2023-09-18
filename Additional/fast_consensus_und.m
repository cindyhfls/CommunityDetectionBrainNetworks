% this code is incomplete because I don't know how to do the triad closing
% part
function ciu = fast_consensus_und(ci,tau,reps,comdet_func,linkidx,m)
% CONSENSUS_UND      Consensus clustering
%
%   CIU = CONSENSUS(D,TAU,REPS) seeks a consensus partition of the 
%   agreement matrix D. The algorithm used here is almost identical to the
%   one introduced in Lancichinetti & Fortunato (2012): The agreement
%   matrix D is thresholded at a level TAU to remove an weak elements. The
%   resulting matrix is then partitions REPS number of times using the
%   Louvain algorithm (in principle, any clustering algorithm that can
%   handle weighted matrixes is a suitable alternative to the Louvain
%   algorithm and can be substituted in its place). This clustering
%   produces a set of partitions from which a new agreement is built. If
%   the partitions have not converged to a single representative partition,
%   the above process repeats itself, starting with the newly built
%   agreement matrix.
%
%   NOTE: In this implementation, the elements of the agreement matrix must
%   be converted into probabilities.
%
%   NOTE: This implementation is slightly different from the original
%   algorithm proposed by Lanchichinetti & Fortunato. In its original
%   version, if the thresholding produces singleton communities, those
%   nodes are reconnected to the network. Here, we leave any singleton
%   communities disconnected.
%
%   Inputs:     
%               TAU,    threshold which controls the resolution of the
%                       reclustering
%               REPS,   number of times that the clustering algorithm is
%                       reapplied
%               
%
%   Outputs:    CIU,    consensus partition
%
%   References: Lancichinetti & Fortunato (2012). Consensus clustering in
%   complex networks. Scientific Reports 2, Article number: 336.
%    Tandon, A., Albeshri, A., Thayananthan, V., Alhalabi, W. & Fortunato, S. Fast consensus clustering in complex networks. Phys. Rev. E 99, 042301 (2019).

%
%   Richard Betzel, Indiana University, 2012
%
%   modified on 3/2014 to include "unique_partitions"
%  modified 9/3/2023 to include fast consensus for big networks, Jiaxin
%  Cindy Tu

n = size(ci,1); 
while 1
    d = sparse(zeros(n));
    dd = agreement(ci)/size(ci,2);
    d(linkidx) = dd(linkidx);
    if mean(nonzeros(dt)<1)<0.02
        break
    end
    % threshold
    dt = d.*(d >= tau).*~eye(n); 
    % triad closing
    mnodeidx = randsample(n,m)';
   
    if nnz(dt) == 0
        ciu = (1:n)';
    else
        ci = zeros(n,reps);
        for iter = 1:reps
            ci(:,iter) = comdet_func(dt);
        end
        ci = relabel_partitions(ci);
        ciu = unique_partitions(ci);
    end
end
ciu = comdet_func(dt);

function cinew = relabel_partitions(ci)
[n,m] = size(ci);
cinew = zeros(n,m);
for i = 1:m
    c = ci(:,i);
    d = zeros(size(c));
    count = 0;
    while sum(d ~= 0) < n
        count = count + 1;
        ind = find(c,1,'first');
        tgt = c(ind);
        rep = c == tgt;
        d(rep) = count;
        c(rep) = 0;
    end
    cinew(:,i) = d;
end

function ciu = unique_partitions(ci)
ci = relabel_partitions(ci);
ciu = [];
count = 0;
c = 1:size(ci,2);
while ~isempty(ci)
    count = count + 1;
    tgt = ci(:,1);
    ciu = [ciu,tgt];                %#ok<AGROW>
    dff = sum(abs(bsxfun(@minus,ci,tgt))) == 0;
    ci(:,dff) = [];
    c(dff) = [];
end