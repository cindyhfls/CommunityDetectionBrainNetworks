% This talks about using silhouette value to choose the optimal number of clusters https://www.mathworks.com/help/stats/clustering.evaluation.silhouetteevaluation.html#bt05vel

function stats=Matrix_metrics_HSB(clrs,rmat0,rth,binarize)
%
% This function calculates the modularity for a set of clustering
% assignments for a given rmatrix. This metric is based on Newman, 2004.
% It is assumed that the clrs per kden/rth are along each colume on clrs.
% Ignore modules set to zero (unclassified).

% Updates:
% 20230530: add in silhouette index, number of communities and number of
% non-singleton communities

%% Set parameters
[Nroi,Nkden]=size(clrs);
stats.modularity=zeros(Nkden,1);
rmat=zeros(size(rmat0),'single');
UDidx=triu(ones(Nroi),1)==1;
domod=1;

%% Calculate modularity
for j=1:Nkden
    
    rmat=rmat0;
    rmat(rmat0<rth(j))=0;                      % Threshold
    if binarize, rmat=single(rmat>0); end      % Binarize
    
    if domod==1
        %         disp('Calculate Modularity');
        clusters=squeeze(clrs(:,j));
        Ugroups=setdiff(unique(clusters),0);
       
        % Initialize
        g=size(Ugroups,1);
        e0=zeros(g,'single');
        
        % Calc strength of edges within and btwn groups
        if isempty(Ugroups)
            stats.modularity = 0;
        else
            for x=1:g
                a=(clusters==Ugroups(x));
                for y=1:g
                    e0(x,y)=nansum(nansum(rmat(a,(clusters==Ugroups(y)))));
                end
            end
            % Calc modularity
            e0=e0/sum(e0(:));
            tr=trace(e0);
            e2=sum((sum(e0)).^2);
            stats.modularity(j)=tr-e2;
        end
    end
    
    %         disp('Calculate Nedges');
    stats.Nedges(j)=sum(rmat(UDidx)~=0);
    %         disp('Calculate kave');
    stats.kave(j)=(stats.Nedges(j)*2)/Nroi;
    %         disp('Calculate Assortativity');
    degrees=sum(rmat)';
    [x,y]=find(triu(rmat,1));
    b=[degrees(x),degrees(y)]-1;% JCT: don't understand why this is -1, it should not affect anything
    c=corrcoef(b);
    stats.A(j)=c(1,2);
    %         disp('Calculate Num Components etc (takes time)');
    M=diag(-sum(rmat))+rmat;    % go go gadget Ncomponents
    [~,D]=eig(M);
    D=abs(diag(D));
    stats.Nc(j)=sum(D(:)<1e-4); % Ncomponents
    D=zeros(Nroi);
    n=1;nPATH=rmat;L=(nPATH~=0);
    while find(L,1), D=D+n.*L;n=n+1;nPATH=nPATH*rmat;L=(nPATH~=0).*(D==0);end
    R=single(D~=0);
    stats.NnBc(j)=max(sum(R))/Nroi; % Nn in biggest component
    stats.Cdns(j)=sum(R(:))/(Nroi*Nroi); % Connectedness
    
 
    M = ones(max(clusters));noneidx = find(unique(clusters)==0);
    M(noneidx,:) = 0; M(:,noneidx) = 0;M = M-diag(diag(M));
    if isempty(noneidx)
        keepnets = true(size(clusters));
    else
        keepnets = clusters~=noneidx;
    end
    D = calc_correlationdist(rmat0);
    % calculate silhouette index
    s = silhouette_coef_mod(clusters(keepnets),D(keepnets,keepnets),M);
    stats.AvgSil(j) = mean(s);
    stats.StdSil(j) = std(s);
    
    % calculate the number of communities
    stats.num_communities(j) = length(setdiff(unique(clusters),0));
    
    % calculate the number of nonsingleton communities
    stats.non_singleton(j) =  count_nonsingleton(clusters(clusters~=0));
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
silhouettes = zeros(n,1);
alternativeid = zeros(n,1);
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

assert(sum(isnan(silhouettes))==0);

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
min_id(ambiguous_ids) = NaN;% update 2023.05.23 set the ambiguous ids to NaN.

if ~isnan(min_id)
    min_id = ids(min_id);
end
end

function [D,R] = calc_correlationdist(zmat,subsample)
% Jiaxin Cindy Tu 2022.11.22
% Calculates the correlation for NxN matrix but excludes the diagonal
% Use subsample to speed up the calculation (for dconns get ~0.99
% correlation for 1/100th subsample
% Unlike the Evan Gordon or Ruby Kong version, the subsampling is not
% random, but more evenly distributed and reproducible by sampling every
% kth vertex

    if ~exist('subsample','var')||isempty(subsample)
        subsample = 1;% don't subsample
    end

    % check that the matrix is the right format
    assert(size(zmat,1)==size(zmat,2));
    assert(length(size(zmat))==2);
    
    zmat(eye(size(zmat))==1) = NaN;% this calculation took too long (~25s)

    R = corr(zmat(1:subsample:end,:),'rows','pairwise');

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

