

%% If loading from saved data
% clear;close all;clc;
% filename = './Infomap/Infomap_BCP_220601.mat'
filename = './Infomap_washu120_low0.001_step0.001_high0.100_xdist20.mat'
load(filename)
% stats = stats{1}; % at sometime in 202205 I changed the Infomap stats results to a cell format so I can save multiple results
params = stats.params

Nroi = size(params.roi,1)

% figdir = fullfile('./Figures',params.IMap_fn);

% load('MNI152nl_on_TT_coord_meshes_32k','MNIl','MNIr'); % adult711B
load('MNI_coord_meshes_32k.mat','MNIl','MNIr');
Anat.CtxL=MNIl;Anat.CtxR=MNIr;
clear MNIl MNIr

%% Label and Color Brain Networks identified by Infomap

% Find the unique networks identified by Infomap
Nets=unique(Cons.SortCons(:));
Nnets=length(Nets);

nameoption =1;
switch nameoption
    case 1 % automatic color
        % % Option 1. % %
        % Auto name 1-#ROI and color based on Jet color lookup table
        AutoName=1;
        CW.cMap = linspecer(Nnets-1);
        CW.Nets=cell(Nnets-1,1);
        for j=1:Nnets    
            if j<10
                CW.Nets{j,1}=['N0',num2str(j)];
            else
                CW.Nets{j,1}=['N',num2str(j)];
            end
        end
        GenOrder=1:max(Cons.SortCons(:));
        
    case 2 % save the network assignment and manually name them
        % % Option 2. % %
        % Name and color networks manually (e.g. for final poster or paper)
        
        % Look at ROIs on cortical surface and take screen shot; press space to advance to next network
        % screen shot networks and save (e.g. to ppt) for labeling
        AutoName=0;
        close all;
        for j=1:Nnets
            Vis_IM_ROI_Module_HSB(Cons.SortCons,stats,Anat,j,Nroi);
            %     print(gcf,['./Figures/',params.IMap_fn,'network',num2str(j)],'-dtiff','-r0');
            img = getframe(gcf);
            %             imwrite(img.cdata,['network',sprintf('%02d',j-1),'.tif'])
            pause;
            close;
        end
        % and then manually update the classification in Util.makeCW
        
    case 3 % name according to template
        template_match_threshold = 0.30; % minimum cutoff for pct_match/sim_mat
        %         template = load('IM_13nets_246.mat');
        %         template = load('IM_246inVol_Talairach_Seitzman2020_14nets(subcortical).mat');
        template =load('IM_Gordon_13nets_333Parcels.mat');
        template.IM.key = sortrows(template.IM.key,1);
        
        nTemplate = max(template.IM.key(:,2));
        nCons = size(Cons.SortCons,2);
        [pct_match,sim_mat] = deal(NaN(nTemplate,Nnets-1,nCons));
        [VIn,MIn] = deal(NaN(1,nCons));
        for iCons = 1:nCons
            tmp = Cons.SortCons(:,iCons);tmp(tmp==0) = find(tmp==0)+1000; % add a large number so 0 is not a single community
            [VIn(iCons), MIn(iCons)] = BCT.partition_distance(template.IM.key(:,2), tmp);
            uniqueMatch = unique(Cons.SortCons(:,iCons));
            uniqueMatch = setdiff(uniqueMatch,0)';
            for i = 1:nTemplate
                for j = uniqueMatch
                    sim_mat(i,j,iCons) = dice(template.IM.key(:,2)==i, Cons.SortCons(:,iCons)==j); % measures how much the overlap between putative network and template
                    pct_match(i,j,iCons) = mean(template.IM.key(Cons.SortCons(:,iCons)==j,2)==i)*100; % measures the percentage division of the nodes in the network belonging to each template
                end
            end
        end
        %         visualize
        figure(991);
        subplot(1,2,1);
        imagesc(nanmean(sim_mat,3));
        xticks(1:Nnets);
        yticks(1:nTemplate);
        yticklabels(template.IM.Nets);
        ytickangle(45);
        xlabel('tentative networks','interpreter','none');
        ylabel(template.IM.name,'interpreter','none');
        colorbar;
        title('dice coefficient');
        subplot(1,2,2);
        imagesc(nanmean(pct_match,3));
        xticks(1:Nnets);
        yticks(1:nTemplate);
        yticklabels(template.IM.Nets);
        ytickangle(45);
        xlabel('tentative networks','interpreter','none');
        ylabel(template.IM.name,'interpreter','none');
        colorbar;
        title('% match to the network');
        
        % print report
        %         [maxv,maxi] = maxk(nanmean(pct_match,3),3);
        %         unclassified = maxv(1,:) < template_match_threshold;
        %         for i = 1:length(maxi)
        %             if unclassified(i)
        %                 fprintf('Network %i: unclassified, %s = %2.0f %%, %s = %2.0f %%, %s = %2.0f %%,\n',i,template.IM.Nets{maxi(1,i)},maxv(1,i),template.IM.Nets{maxi(2,i)},maxv(2,i),template.IM.Nets{maxi(3,i)},maxv(3,i));
        %             else
        %                 fprintf('Network %i: %s, %s = %2.0f %%, %s = %2.0f %%, %s = %2.0f %%,\n',i,template.IM.Nets{maxi(1,i)},template.IM.Nets{maxi(1,i)},maxv(1,i),template.IM.Nets{maxi(2,i)},maxv(2,i),template.IM.Nets{maxi(3,i)},maxv(3,i));
        %             end
        %         end
        
        [maxv,maxi] = maxk(nanmean(sim_mat,3),3);
        unclassified = maxv(1,:) < template_match_threshold;
        for i = 1:length(maxi)
            if unclassified(i)
                fprintf('Network %i: unclassified, %s = %1.2f , %s = %1.2f , %s = %1.2f ,\n',i,template.IM.Nets{maxi(1,i)},maxv(1,i),template.IM.Nets{maxi(2,i)},maxv(2,i),template.IM.Nets{maxi(3,i)},maxv(3,i));
            else
                fprintf('Network %i: %s, %s = %1.2f , %s = %1.2f , %s = %1.2f ,\n',i,template.IM.Nets{maxi(1,i)},template.IM.Nets{maxi(1,i)},maxv(1,i),template.IM.Nets{maxi(2,i)},maxv(2,i),template.IM.Nets{maxi(3,i)},maxv(3,i));
            end
        end
        
        % and then manually update the classification in Util.makeCW
    case 4 % already have the names saved in another file (after 2) and loading that
        % Generate Consensus-wide naming/coloring convention structure
        [CW,GenOrder] = Util.makeCW(params,Nnets); % I stored my manual assignment in this function
        % N.B. the network names are in /DIAN folder hmmm
end



%% Re-Order Networks (vis,DMN,Mot,DAN,FPC,...)

% % Option 1. Specify the order your networks will appear in figures
% If you named and colored your networks, you may want them to appear in an
% order mirrioring published papers (e.g. Eggebrecht et al., 2017 or Power et al., 2011)
% orig_order=[1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20];
% GenOrder = [1,8,13,2,4, 9,10,18,5,12, 16,7,3,6,17, 11,19,14,15,20]';
% GenOrder = [3,7,11,2,13,12,1,6,4,10,8,9,5];
% GenOrder = [4,2,1,6,3,5,7:13];

% % Option 2. Leave network order as default order (typically default output is largest to smallest)
% GenOrder=1:max(Cons.SortCons(:));

% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=Cons.SortCons;foo(foo==0)=NaN;
Cons.SortConsRO=Cons.SortCons;
for j=1:length(GenOrder),Cons.SortConsRO(foo==GenOrder(j))=j;end
foo=Cons.SortConsRO;


%% Visualize Networks-on-brain and Consensus Edge Density Matrix
Explore_ROI_kden_HSB(foo,CWro.cMap,Anat,params.roi,Cons.epochs.mean_kden);

%% Generate Infomap (IM) Structure for Viable Edge Density Ranges
% Viable = edge densities in which connectivity >80% (see figure output
% from Org_Cons_Org_Imap_Matrix)
% IM structures are used during Enrichment to organize ROI into networks
% This set of codes visualizes all possible IM options to choose from

% load FC matrix
if strcmp(stats.params.format,'mat')
    stats.MuMat = smartload(stats.params.zmatfile); %(parcel-wise data in .mat)
elseif strcmp(stats.params.format,'cifti')
    tmp = ft_read_cifti_mod(stats.params.zmatfile);stats.MuMat = tmp.data; %(vertex-wise datain cifti format)
end

toIM=[1:size(Cons.SortCons,2)];
for j=toIM % Auto out of IM for each Cons model
    
    %     if (Cons.stats.NnBc(j)>0.9) && (Cons.stats.kave(j)>log(Nroi))
    % remove number of nodes in largest component <= 90%? and mean degree
    % <log(Nroi)? degree is calculated with number of non-zero weight connections
    
    
    % Turn into function? inputs: Cons, CWro, IM name, stats
    
    cMap=CWro.cMap;
    Nets=CWro.Nets;
    
    % Add a way to fix USp?
    temp=Fix_US_Grouping_HSB(Cons.SortConsRO,j); %This code attempts to assign unspecified ROIs to networks, when possible
    %         temp(string(Nets)=='None'|string(Nets)=='Usp')=0;
    % temp=squeeze(Cons.SortConsRO(:,j));
    
    % USp networks with less than 5
    NnetsI=unique(temp);
    for nn=1:length(NnetsI)
        if sum(temp==NnetsI(nn))<5,temp(temp==NnetsI(nn))=0;end
    end
    
    if any(temp==0)
        temp(temp==0)=size(cMap,1)+1;
        cMap=cat(1,cMap,[0.25,0.25,0.25]);% gray for USp
        %             cMap = cat(1,cMap,[1,1,0.8]); % a very light yellow for USp
        Nets=cat(1,Nets,'None');
    end
    keep=unique(temp)';
    
    % Put together IM structure
    clear IM
    %     IM.name = ['IM_',params.IMap_fn,'_Consesus_model_winnertakesall'];
    IM.name=['IM_',params.IMap_fn,'_Consesus_model_',num2str(j)];
    IM.cMap=cMap(keep,:);
    IM.Nets=Nets(keep);
    IM.ROIxyz=params.roi;
    IM.key=[[1:Nroi]',zeros(Nroi,1)];
    [IM.key(:,2),IM.order]=sort(IM_Remove_Naming_Gaps_HSB(temp));
    IM.ROIxyz=IM.ROIxyz(IM.order,:);
    IM=Org_IM_DVLR_HSB(IM);    
    
    % Visualize
    figure('Color','k','Units','Normalized','Position',[0.35,0.30,0.35,0.61]);
    subplot(4,1,[1:3])
    Matrix_Org3_HSB(stats.MuMat(IM.order,IM.order),IM.key,10,[-0.6,0.6],IM.cMap,0);
    title([{[strrep(IM.name,'_',' ')]};{['kden=',...
        num2str(stats.kdenth(Cons.epochs.kden_i(j))*100,'%0.2f'),...
        '-',num2str(stats.kdenth(Cons.epochs.kden_f(j))*100,'%0.2f'),'%, ',...
        num2str(length(IM.Nets)),' Networks']}],'Color','w')
    c=colorbar;c.Ticks=[-0.6,0,0.6];c.TickLabels={'-0.6','0','0.6'};
    set(c,'Color','w')
    subplot(4,1,4)
    histogram(IM.key(:,2));title(strrep(IM.name,'_',' '))
    set(gca,'XTick',[1:max(IM.key(:,2))],'XTickLabel',...
        IM.Nets,'Color','k',...
        'XColor','w','YColor','w');
    xtickangle(45);
    ylabel('Nrois','Color','w')
    xlim([0,max(IM.key(:,2)+1)]);
    set(gcf,'InvertHardCopy','off');
%     print(gcf,fullfile(figdir,sprintf('%s_heatmap',IM.name)),'-dtiff')
    
    % if function, end here
    
    % Save the IM file
%     save(['./Results/',IM.name],'IM');
    
    %     end
end

%% Visualize Mean and Variance for IM model (double click IM model to load)
% figure('Color','w','Units','Normalized','Position',[0.4,0.4,0.5,0.4]);
% subplot(1,2,1);imagesc(rmatAve(IM.order,IM.order),[-0.6,0.6]);
% axis square;title('Mean');colormap(gca,'jet');colorbar;
% subplot(1,2,2);imagesc(rmatVar(IM.order,IM.order));
% axis square;title('Variance');colormap(gca,'parula');colorbar;

%% Visualize IM Model on Brain with Network Names and Colors
params.radius = 4;
Anat.alpha = 1;
IM.cMap(end,:) = [0.3,0,0.6];
Matrix_Org3_HSB(stats.MuMat(IM.order,IM.order),...
    IM.key,10,[-0.3,0.3],IM.cMap,1); % mean
% Anat.ctx='std';Plot.View_ROI_Modules(IM,Anat,IM.ROIxyz,params);
% figure;histogram(IM.key(:,2));title(strrep(IM.name,'_',' '))
% set(gca,'XTick',[1:max(IM.key(:,2))],'XTickLabel',IM.Nets);ylabel('Nrois')
%%
to_load = 'IM_Infomap_BCP_gp5_220304_Consesus_model_4';IM1 = load(to_load);
left_lab = IM1.IM.Nets;
to_load = 'IM_Infomap_BCP_gp4_220307_Consesus_model_4';IM2 = load(to_load);
right_lab = IM2.IM.Nets;

D = NaN(length(left_lab),length(right_lab));
for i = 1:length(left_lab)
    for j = 1:length(right_lab)
        D(i,j) = sum(IM2.IM.key(IM1.IM.key(:,2)==i,2)==j);
    end
end

h = Plot.alluvialflow(D, left_lab, right_lab, 'gp5(Consensus4) - gp4 (Consensus 4)',IM1.IM.cMap)

