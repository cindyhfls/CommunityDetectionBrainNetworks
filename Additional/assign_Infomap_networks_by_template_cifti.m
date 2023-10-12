function [CW,GenOrder,MIn] = assign_Infomap_networks_by_template_cifti(Cons,template_cifti,template_match_threshold,template_match_method)
% this function takes a template and make that match to the infomap result
% assumes the input cifti is a dlabel file with labels
Nets=setdiff(unique(Cons.SortCons(:)),0);
Nnets=length(Nets);
NetNames = {template_cifti.diminfo{2}.maps.table(2:end).name}';
NetcMap =horzcat(template_cifti.diminfo{2}.maps.table(2:end).rgba)';NetcMap = NetcMap(:,1:3);
if ~exist('template_match_threshold','var')||isempty(template_match_threshold)
    template_match_threshold = 0;
end
if ~exist('template_match_method','var')||isempty(template_match_method)
    template_match_method = 'dice';
end
%% Find the MI and overlap for each network
        templateKey =template_cifti.cdata;

        nTemplate = max(templateKey);
        nCons = size(Cons.SortCons,2);
        [pct_match,sim_mat] = deal(NaN(nTemplate,Nnets-1,nCons));
        [VIn,MIn] = deal(NaN(1,nCons));
        for iCons = 1:nCons
            tmp = Cons.SortCons(:,iCons);
            %tmp(tmp==0) = find(tmp==0)+1000; % add a large number so 0 is not a single community
            idx = tmp~=0 &templateKey~=0;
            [VIn(iCons), MIn(iCons)] = partition_distance(templateKey(idx), tmp(idx));
            uniqueMatch = unique(Cons.SortCons(:,iCons));
            uniqueMatch = setdiff(uniqueMatch,0)';
            for i = 1:nTemplate
                for j = uniqueMatch
                    sim_mat(i,j,iCons) = dice(templateKey==i, Cons.SortCons(:,iCons)==j); % measures how much the overlap between putative network and template
                    pct_match(i,j,iCons) = mean(templateKey(Cons.SortCons(:,iCons)==j)==i)*100; % measures the percentage division of the nodes in the network belonging to each template
                end
            end
        end
        %% visualize
        figure(991);
        switch template_match_method
            case 'dice'
                imagesc(nanmax(sim_mat,[],3));
%                 imagesc(nanmean(sim_mat,3));
                title('dice coefficient');
            case 'percentage'
                imagesc(nanmax(pct_match,[],3));
%                 imagesc(nanmean(pct_match,3));
                title('% composition of network');
        end
        xticks(1:Nnets);
        yticks(1:nTemplate);
        yticklabels(NetNames)
        ytickangle(45);
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
        CW.Nets =NetNames(idx);
        CW.Nets(maxv(1,:)<template_match_threshold) = repelem({'None'},sum(maxv(1,:)<template_match_threshold),1);

        CW.cMap = NetcMap(idx,:);
        CW.cMap(maxv(1,:)<template_match_threshold,:) = 0.5;
        
        uniqueNets = setdiff(unique(CW.Nets),{'None','USp'});
        
        for i = 1:length(uniqueNets)
            matchnet = find(string(CW.Nets)==uniqueNets{i});
            if length(matchnet)>1
                for k = 2:length(matchnet)
                    CW.Nets{matchnet(k)}=[CW.Nets{matchnet(k)},num2str(k)];
                    CW.cMap(matchnet(k),:) = change_rgb_color(CW.cMap(matchnet(k-1),:));
                end
            end
        end
        
        idxwithnone = idx;
        idxwithnone(any(string(CW.Nets)==["None","USp"],2))= Inf;
        [~,GenOrder] = sort(idxwithnone);
   
end