function [CW,GenOrder] = assign_networks_by_template_cifti(Clust,template_cifti,template_match_threshold,template_match_method)
% this function takes a template and make that match to the infomap result
% assumes the input cifti is a dlabel file with labels
Nets=setdiff(unique(Clust(:)),0);
Nnets=length(Nets);
NetNames = {template_cifti.diminfo{2}.maps.table(2:end).name}';
NetcMap =horzcat(template_cifti.diminfo{2}.maps.table(2:end).rgba)';NetcMap = NetcMap(:,1:3);

if ~exist('template_match_threshold','var')||isempty(template_match_threshold)
    template_match_threshold = 0;
end
if ~exist('template_match_method','var')||isempty(template_match_method)
    template_match_method = 'dice';
end
if strcmp(template_match_method,'percentage')
    template_match_threshold = template_match_threshold*100;
end
%% Use the most representative network
G1=setdiff(unique(Clust(:)),0);
% for j=1:length(G1)
%     tmp =sum(Clust==G1(j,1));tmp(tmp==0) = NaN;
%     [G1(j,2)]=max(tmp); % Adam used max when making manual judgement of network names
%     G1(j,3) = find(tmp==G1(j,2),1);
% end
% repnets = cell2mat(arrayfun(@(ii)Clust(:,G1(ii,3))==G1(ii,1),G1(:,1),'UniformOutput',false)');

%% Find the  overlap for each network
Nnets = length(G1);
templateKey =template_cifti.cdata;
nTemplate = max(templateKey);
nlevels = size(Clust,2);
[pct_match,sim_mat] = deal(NaN(nTemplate,Nnets,nlevels));

for j = 1:Nnets
    j
    for i = 1:nTemplate
        for k = 1:nlevels
           sim_mat(i,j,k) = dice((templateKey==i),Clust(:,k)==j);
           pct_match(i,j,k) = mean(templateKey(Clust(:,k)==j)==i)*100;
        end
%         sim_mat(i,j) = dice((templateKey==i),repnets(:,j));
%         pct_match(i,j) = mean(templateKey(repnets(:,j))==i)*100;
    end
end
%% visualize
figure(991);
switch template_match_method
    case 'dice'
        imagesc(nanmax(sim_mat,[],3));
        title('dice coefficient');
    case 'percentage'
        imagesc(nanmax(pct_match,[],3));
        title('% composition of network');
end
xticks(1:Nnets);
yticks(1:nTemplate);
yticklabels(NetNames)
ytickangle(45);
xtickangle(90);
xlabel('tentative networks','interpreter','none');
colorbar;
set(gca,'TickLabelInterpreter','none');
%% Display Output
switch template_match_method
    case 'dice'
        %       print report for dice
        disp('% dice similarity to network');
        [maxv,maxi] = maxk(nanmax(sim_mat,[],3),3);
        unclassified = maxv(1,:) < template_match_threshold;
        for i = 1:length(maxi)
            if unclassified(i)
                fprintf('Network %i: unclassified, %s = %1.2f , %s = %1.2f , %s = %1.2f ,\n',i,NetNames{maxi(1,i)},maxv(1,i),NetNames{maxi(2,i)},maxv(2,i),NetNames{maxi(3,i)},maxv(3,i));
            else
                fprintf('Network %i: %s, %s = %1.2f , %s = %1.2f , %s = %1.2f ,\n',i,NetNames{maxi(1,i)},NetNames{maxi(1,i)},maxv(1,i),NetNames{maxi(2,i)},maxv(2,i),NetNames{maxi(3,i)},maxv(3,i));
            end
        end
    case 'percentage'
        disp('% composition of network');
        [maxv,maxi] = maxk(nanmax(pct_match,[],3),3);
        unclassified = maxv(1,:) < template_match_threshold;
        for i = 1:length(maxi)
            if unclassified(i)
                fprintf('Network %i: unclassified, %s = %2.0f %%, %s = %2.0f %%, %s = %2.0f %%,\n',i,NetNames{maxi(1,i)},maxv(1,i),NetNames{maxi(2,i)},maxv(2,i),NetNames{maxi(3,i)},maxv(3,i));
            else
                fprintf('Network %i: %s, %s = %2.0f %%, %s = %2.0f %%, %s = %2.0f %%,\n',i,NetNames{maxi(1,i)},NetNames{maxi(1,i)},maxv(1,i),NetNames{maxi(2,i)},maxv(2,i),NetNames{maxi(3,i)},maxv(3,i));
            end
        end
end
%% Save to output
idx = maxi(1,:);
vals = maxv(1,:);
CW.Nets =NetNames(idx);
CW.Nets(maxv(1,:)<template_match_threshold) = repelem({'USp'},sum(maxv(1,:)<template_match_threshold),1);

CW.cMap = NetcMap(idx,:);

uniqueNets = setdiff(unique(CW.Nets),{'None','USp'});
validNets = any(string(CW.Nets)==string(uniqueNets)',2);
n_new_color = sum(validNets)-length(uniqueNets);
n_max = 150; % maximum possible distinguishable colors
assert(n_max>n_new_color,'too many unique clusters, have you really matched across columns and removed the singleton ones?'); % error out if too many colors
% new_color = distinguishable_colors(n_max-n_new_color,unique(CW.cMap,'rows')); % we should just set it to new_color = distinguishable_colors(n_max)
new_color = distinguishable_colors(n_max);
mind_toexisting= min(pdist2(rgb2lab(NetcMap),rgb2lab(new_color),'Euclidean'));% remove the colors in the template from new_color
new_color(mind_toexisting<20,:) = [];

for i = 1:length(uniqueNets)
    matchnet = find(string(CW.Nets)==uniqueNets{i});
    [~,sortmatch] = sort(vals(matchnet),'descend'); % update 2023/12/26 instead of giving the first network the core color, giving it to the maximum match one instead
    matchnet = matchnet(sortmatch);
    n_colors = length(matchnet);
    if  n_colors >1
        core_clr = CW.cMap(matchnet(1),:);
        for k = 2:length(matchnet)
            % for each new color assignment, first exclude all
            % possibilities that are too similar to the existing colors
            mind_toexisting= min(pdist2(rgb2lab(unique(CW.cMap,'rows')),rgb2lab(new_color),'Euclidean'));
            new_color(mind_toexisting<20,:) = [];
            
            % then find the color that is the most similar to the core
            % color and assign to the second matching network
            d = pdist2(rgb2lab(core_clr),rgb2lab(new_color),'Euclidean');
            if isempty(d)
                error('too many colors to assign')
            end
            [~,minid] = min(d);
            
            CW.Nets{matchnet(k)}=[CW.Nets{matchnet(k)},'_',num2str(k)];
            CW.cMap(matchnet(k),:) = new_color(minid,:);%change_rgb_color(CW.cMap(matchnet(k-1),:),CW.cMap);
            
            % finally remove that color from being considered as a
            % possibility
            new_color(minid,:) = [];
        end
    end
end

% now match for the unspecified networks
matchnet = find(string(CW.Nets)=='USp');
n_USp = length(matchnet);

for k= 1:n_USp
    d = pdist2(rgb2lab(new_color),rgb2lab(unique(CW.cMap,'rows')),'Euclidean');
    mind = min(d,[],2);
    [~,ids] = max(mind); % find n_USp new colors that are most different from the existing ones
    CW.Nets{matchnet(k)}=[CW.Nets{matchnet(k)},'_',num2str(k)];
    CW.cMap(matchnet(k),:) = new_color(ids,:);%change_rgb_color(CW.cMap(matchnet(k-1),:),CW.cMap);
end

idxwithnone = idx;
idxwithnone(contains(string(CW.Nets),["None","USp"]))= Inf;
[~,GenOrder] = sort(idxwithnone);

end