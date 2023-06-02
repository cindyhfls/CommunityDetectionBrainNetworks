function [Y E] = cmdscale_matHSB(mat,groups,distmet,varargin)
% This function performs multi-dimensional scaling on input matrices and
% displays a 2-dimensional plot of the relative positions of the input 
% matrices colored by group. MDS is performed using euclidean distance.
% Inputs:
% mat - node x node x subject array of matrices from all subjects
% groups - identifies which subject belongs to which group
% varargin - group names can be specified by a third input, e.g. 
% {'group1';'group2'} and will be displayed in a legend.
% Outputs:
% Y - configuration matrix
% E - eigenvalues of Y*Y' 
% 
% TOL, 09/14
numnodes = size(mat,1);
numgroups = length(unique(groups));
colors = [1,0,0;0,0,1;0,0.5,0;1,1,0;0,0,0;...%r,b,g,y,k
    0.8,0.6,0;0.75,0,0.75;0,1,1;0.4,0.2,0;0.2,1,0.2];%o,m,c,b,lg
if nargin > 3
    groupnames = varargin{1};
end

mask = ones(numnodes);
mask = triu(mask,1);
for s = 1:size(mat,3)
        temp = mat(:,:,s);
        mat_col(:,s) = temp(logical(mask));
end

% Multi-dimensional scaling
D = pdist(double(mat_col'),distmet);
[Y E] = cmdscale(D);

% Display result
figure('Color','white','Position',[100,100,800,750])    
hold
% individuals
for c = 1:numgroups
    inds = groups==c;
    p{c}=plot(Y(inds,1),Y(inds,2),'o',...
        'Color',colors(c,:),'MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(c,:));
end
% Group ave
for c = 1:numgroups
    inds = groups==c;
    
    plot(mean(Y(inds,1)),mean(Y(inds,2)),'^','Color',colors(c,:),...
        'MarkerSize',16,...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(c,:))
    % add SE cross
    plot([mean(Y(inds,1))-std(Y(inds,1))/sqrt(sum(inds)),...
        mean(Y(inds,1))+std(Y(inds,1))/sqrt(sum(inds))],...
        [mean(Y(inds,2)),mean(Y(inds,2))],...
        'k','LineWidth',2)%'Color',colors(c,:),
    plot([mean(Y(inds,1)),mean(Y(inds,1))],...
        [mean(Y(inds,2))-std(Y(inds,2))/sqrt(sum(inds)),...
        mean(Y(inds,2))+std(Y(inds,2))/sqrt(sum(inds))],...
        'k','LineWidth',2)%,'Color',colors(c,:)
end
set(gca,'FontWeight','bold','FontSize',14)
if nargin > 2
    legend([p{:}],groupnames,'FontWeight','bold','FontSize',14,...
        'Location','best')
end
xlabel('MDS dimension 1');ylabel('MDS dimension 2')
axis square