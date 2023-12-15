
function stats=Matrix_metrics_HSB(clrs,rmat0)
%
% This function calculates some basic summary metrics of the communities
% Input: 
% clrs - cluster assignment in nxp where n is the number of nodes
% and p is the levels (e.g. different K for k-means, different layers for
% multilayer networks). p can also be repeats, p can also 1.
% rmat0 - the matrix in nxn

% Updates:
% 20230530: add in silhouette index, number of communities and number of
% non-singleton communities

% 20230914 change input structure

% 20231127 remove the legacy code way of calculating Newman modularity on
% thresholded matrix and the connectedness measures etc.

% This post talks about using silhouette value to choose the optimal number of clusters https://www.mathworks.com/help/stats/clustering.evaluation.silhouetteevaluation.html#bt05vel

%% make sure rmat0 has no diagonal element
for ii = 1:length(rmat0)
    rmat0(ii,ii) = 0;
end

%% Set parameters
[Nroi,Nlevels]=size(clrs);

for j = 1:Nlevels
    clusters=squeeze(clrs(:,j));
   %% Basic stats
    % number of communities and nonsingleton communities
    stats.non_singleton(j) =  count_nonsingleton(clusters(clusters~=0));
    stats.num_communities(j) = length(setdiff(unique(clusters),0));
    
   %% Cluster validity
    M = ones(max(clusters)); M = M-diag(diag(M)); % comparing every community to every other community, excluding it self
    keepnets = clusters~=0; % remove the zero assignments (unassigned community label)
    D = calc_correlationdist(rmat0); % same as pdist(rmat0,'corr') except for ignoring the diagonal
    
    % calculate silhouette index (can also use the matlab silhouette but that may count the diagonals)
    s = silhouette_coef_mod(clusters(keepnets),D(keepnets,keepnets),M); % use my own silhouette instead of MATLAB because I wanted to calculate correlation distance without the diagonal and also MATLAB R2020b version seems to get the wrong numbers
    stats.AvgSil(j) = nanmean(s);
    stats.StdSil(j) = nanstd(s);
    
    % Calculate Davies-Bouldin index
    tmp = evalclusters(double(rmat0(keepnets,keepnets)),clusters(keepnets),'DaviesBouldin');
    stats.DB(j) = tmp.CriterionValues;
    
    % Calculate Calinski-Harabasz index
    tmp = evalclusters(double(rmat0(keepnets,keepnets)),clusters(keepnets),'CalinskiHarabasz');
    stats.CH(j) = tmp.CriterionValues;
    
    %% Modularity
    % use the Rubinov & Sporns 2011 modularity on signed networks
    stats.Qsigned(j) = modularity_signed(rmat0,clusters);
    stats.QNG(j) =modularity_newman(rmat0,clusters);
    %% System segregation
    % Reference: Chan et al. 2014 PNAS Decreased segregation of brain
    % systems across the healthy adult lifespan
    within= rmat0(bsxfun(@eq,clusters,clusters.'));
    Zw = mean(nonzeros(within));
    between = rmat0(bsxfun(@ne,clusters,clusters.'));
    Zb = mean(nonzeros(between));
    stats.SysSeg(j) = (Zw-Zb)/Zw;    
     
    % calculate silhouette index with the correlation matrix itself already
    % as a distance measure (temporal similarity instead of connectivity
    % spatial similarity)
    D = 1-rmat0; for ii = 1:length(D),D(ii,ii)=0;end
    [avg_within,min_avg_between] = deal(NaN(Nroi,1));
    for iclu=unique(clusters)'
        avg_within(clusters==iclu,1) = sum(D(clusters==iclu,clusters==iclu),2)/sum(clusters==iclu)-1;
        for jclu = setdiff(unique(clusters)',iclu)
            tmp = sum(D(clusters==iclu,clusters==jclu),2)/sum(clusters==jclu);
        end
        min_avg_between(clusters==iclu,1) = min(tmp);
    end
    stats.SilTemporal(j)=mean((min_avg_between-avg_within)/max(min_avg_between-avg_within));
end


end

function [ silhouettes,alternativeid ] = silhouette_coef_mod( parcels, D, neigh )
%SILHOUETTE_COEF Silhouette coefficient of a parcellation.
%   For each vertex in a cortical surface, SILHOUETTE_COEF compares the
%   within-parcel dissimilarity defined as the average distance to all
%   vertices in the same parcel, to the inter-parcel dissimilarity
%   obtained from those assigned to other parcels.
%
%   Let Pi be the parcel to which vertex i is assigned and ai and bi be the
%   average distance from i to the vertices in Pi and to the vertices in
%   other parcels adjacent to Pi. Silhouette coefficient (SI) for i is then
%   computed as:
%
%   SI(i) = (bi - ai)/max(ai,bi)
%
%   This guarantees SI values within [-1, +1], as long as a distance
%   measure is used.
%   Adapted from silhouette_coef.m from Arslan et al. 2018 NeuroImage.

n = length(D);
num = max(parcels);
silhouettes = NaN(n,1);
alternativeid = NaN(n,1);
for i = 1 : num
    in_members = parcels == i;
    nk = sum(in_members);
    
    if nk < 2
        continue; % Singleton parcel detected. Can happen with N-Cuts.
    end
    
    dists = D(in_members,in_members);
    dists(logical(eye(length(dists)))) = 0;
    dists_in = sum(dists,2)/(length(dists)-1);
    
    ids = find(neigh(:,i));
    [dists_out,min_id] = find_dist_to_neighs( ids, parcels, D, in_members );

    silhouettes(in_members) = (dists_out - dists_in) ./ ...
        max([dists_in dists_out],[],2);
    
    
    alternativeid(in_members) = min_id;
end


end


% Compute distance to/from vertices in adjacent parcels
function [ min_votes_out,min_id] = find_dist_to_neighs(ids, parcels, D, in_members )

% Jiaxin Cindy Tu 20221121 the Arslan method use mean across any neighbors
% which is kind of weird! should be the minimum

votes_out = NaN(sum(in_members),length(ids));
for j = 1 : length(ids)
    out_members = parcels == ids(j);
    nk = sum(out_members);
    corrs = D(in_members,out_members);
    votes_out(:,j) = sum(corrs,2)/nk;
end
[min_votes_out,min_id] = min(votes_out,[],2);

ambiguous_ids = sum(min_votes_out==votes_out,2)~=1;
min_id(ambiguous_ids) = NaN; % update 2023.05.23 set the ambiguous ids to NaN.

if ~isnan(min_id)
    min_id = ids(min_id);
end
end

function [D,R] = calc_correlationdist(zmat)
% Jiaxin Cindy Tu 2022.11.22
    % check that the matrix is the right format
    assert(size(zmat,1)==size(zmat,2));
    assert(length(size(zmat))==2);
    
    for i = 1:length(zmat)
        zmat(i,i) = NaN;
    end
    
    R = corr(zmat,'rows','pairwise');
    D = 1-R;
    
end

function nonsingleton = count_nonsingleton(x)
    xs = sort(x);
    counter = 0;
    for i = 2:length(xs)-1
        counter=counter+((xs(i)~=xs(i+1))&&(xs(i)~=xs(i-1)));
    end
    nonsingleton = length(unique(xs))-counter;
end

function Q = modularity_signed(W,M, gamma)
% Rubinov & Sporns 2011
% adapted from community_louvain.m in BCT toolbox
%   Inputs:
%       W,
%           directed/undirected weighted/binary connection matrix with
%           positive and possibly negative weights.
%       M,
%           community affiliation vector
%       gamma,
%           resolution parameter (optional)
%               gamma>1,        detects smaller modules
%               0<=gamma<1,     detects larger modules
%               gamma=1,        classic modularity (default)
if ~exist('gamma','var')||isempty(gamma)
    gamma=1;        % classic modularity (default)
end

W=double(W);                                % convert to double format

W0 = W.*(W>0);                          %positive weights matrix
s0 = sum(sum(W0));                      %weight of positive links
B0 = W0-gamma*(sum(W0,2)*sum(W0,1))/s0; %positive modularity

W1 =-W.*(W<0);                          %negative weights matrix
s1 = sum(sum(W1));                      %weight of negative links
if s1                                   %negative modularity
    B1 = W1-gamma*(sum(W1,2)*sum(W1,1))/s1;
else
    B1 = 0;
end

B = B0/s0 - B1/(s0+s1);


B = (B+B.')/2;                                          % symmetrize modularity matrix

Q = sum(B(bsxfun(@eq,M,M.')));                        % compute modularity
end

function Q = modularity_newman(A,M)
deg = sum(A,2);
twom = sum(deg);
P = (deg*deg')/twom;
B = A - P;
B = (B+B.')/2;                                          % symmetrize modularity matrix
Q = sum(B(bsxfun(@eq,M,M.')))/twom;                        % compute modularity
end