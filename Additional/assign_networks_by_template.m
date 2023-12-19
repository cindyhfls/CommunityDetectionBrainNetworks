function [CW,GenOrder] = assign_networks_by_template(Clust,template,template_match_threshold,template_match_method)
% this function takes a template and make that match to the infomap result
Nets=setdiff(unique(Clust(:)),0);
Nnets=length(Nets);
NetNames = template.IM.Nets;
NetcMap =template.IM.cMap;
if ~exist('template_match_threshold','var')||isempty(template_match_threshold)
    template_match_threshold = 0;
end
if ~exist('template_match_method','var')||isempty(template_match_method)
    template_match_method = 'dice';
end
%% Merge the ones with high overlap (sometimes they exist because breaks in the assignment?)
G1=setdiff(unique(Clust(:)),0);
for j=1:length(G1)
    [G1(j,2),G1(j,3)]=max(sum(Clust==G1(j,1)));
end
repnets = cell2mat(arrayfun(@(ii)Clust(:,G1(ii,3))==G1(ii,1),G1(:,1),'UniformOutput',false)');
%% Find the  overlap for each network
Nnets = size(repnets,2)
[~,sortid] = sort(template.IM.order);
templateKey =template.IM.key(sortid,2);
nTemplate = max(templateKey);
[pct_match,sim_mat] = deal(NaN(nTemplate,Nnets));

for j = 1:Nnets
    for i = 1:nTemplate
        sim_mat(i,j) = dice((templateKey==i),repnets(:,j));
        pct_match(i,j) = mean(templateKey(repnets(:,j))==i)*100;
    end
end

%% visualize
figure(991);
switch template_match_method
    case 'dice'
        imagesc(nanmax(sim_mat,[],3));
        xticks(1:Nnets);
        yticks(1:nTemplate);
        yticklabels(NetNames);
        ytickangle(45);
        xlabel('tentative networks','interpreter','none');
        ylabel(template.IM.name,'interpreter','none');
        colorbar;
        title('dice coefficient');
    case 'percentage'
        imagesc(nanmax(pct_match,[],3));
        xticks(1:Nnets);
        yticks(1:nTemplate);
        yticklabels(NetNames);
        ytickangle(45);
        xlabel('tentative networks','interpreter','none');
        ylabel(template.IM.name,'interpreter','none');
        colorbar;
        title('% composition of network');
end
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
CW.Nets =NetNames(idx);
CW.Nets(maxv(1,:)<template_match_threshold) = repelem({'None'},sum(maxv(1,:)<template_match_threshold),1);

CW.cMap = NetcMap(idx,:);
CW.cMap(maxv(1,:)<template_match_threshold,:) = 0.5;

% uniqueNets = setdiff(unique(CW.Nets),{'None','USp'});
uniqueNets=unique(CW.Nets);
validNets = any(string(CW.Nets)==string(uniqueNets)',2);
n_new_color = sum(validNets)-length(uniqueNets);

new_color = distinguishable_colors(n_new_color,unique(CW.cMap,'rows'));

for i = 1:length(uniqueNets)
    matchnet = find(string(CW.Nets)==uniqueNets{i});
    n_colors = length(matchnet);
    if  n_colors >1
        core_clr = CW.cMap(matchnet(1),:);
        d = pdist2(rgb2lab(core_clr),rgb2lab(new_color),'Euclidean');
        [~,minid] = mink(d,n_colors-1); % take the one that is relatively closer in perception
        while minval<20% prevent being too similar
            new_color(minid,:) = [];
            d(minid) = [];
            [minval,minid] = mink(d,n_colors-1); % take the one that is relatively closer in perception
        end
        for k = 2:length(matchnet)
            CW.Nets{matchnet(k)}=[CW.Nets{matchnet(k)},'_',num2str(k)];
            CW.cMap(matchnet(k),:) = new_color(minid(k-1),:);%change_rgb_color(CW.cMap(matchnet(k-1),:),CW.cMap);
        end
        new_color(minid,:) = [];
    end
end

idxwithnone = idx;
idxwithnone(contains(string(CW.Nets),["None","USp"]))= Inf;
[~,GenOrder] = sort(idxwithnone);

end