function [CWro,SortClusRO] = assign_network_colors(SortClus,nameoption,templatepath,parcelpath)
%% Set default template
if  nameoption==3
    if~exist('templatepath','var')||isempty(templatepath)
        templatepath = 'IM_Gordon_2014_333_Parcels.mat';
    end
end

%% Label and Color Brain Networks identified by Infomap

% Find the unique networks identified by Infomap
Nets=setdiff(unique(SortClus(:)),0);
Nnets=length(Nets);

% nameoption = 3
switch nameoption
    case 1 % automatic color
        % % Option 1. % %
%         CW.cMap = linspecer(Nnets-1);
        CW.cMap = distinguishable_colors(Nnets);
        CW.Nets=cell(Nnets,1);
        for j=1:Nnets    
            CW.Nets{j,1} = sprintf('%02d',j);
        end
        GenOrder=1:max(SortClus(:));
        
    case 2 % save the network assignment and manually name them - this needs some work as I changed the input structure
        load('MNI_coord_meshes_32k.mat')
        Anat.CtxL = MNIl;Anat.CtxR = MNIr;
        clear MNIl MNIr
        % % Option 2. % %
        % Name and color networks manually (e.g. for final poster or paper)
        
        % Look at ROIs on cortical surface and take screen shot; press space to advance to next network
        % screen shot networks and save (e.g. to ppt) for labeling
%         close all;
%         for j=1:Nnets
%             Vis_IM_ROI_Module_HSB(Cons.SortClus,stats,Anat,j,Nroi);
%             pause;
%             close;
%         end
        % and then manually update the classification inmakeCW
        
    case 3 % name according to a given template
        if isstruct(templatepath) % e.g. colortemplate.IM
            colortemplate = templatepath;
            [CW,GenOrder] = assign_networks_by_template(SortClus,colortemplate,0.1,'dice');%'dice'
        elseif contains(templatepath,'.mat') % e.g. load IM structure
            colortemplate =load(templatepath);
            [CW,GenOrder] = assign_networks_by_template(SortClus,colortemplate,0.1,'dice');%'dice'
        elseif contains(templatepath,'.nii') % e.g. load cifti file
            if isnumeric(parcelpath)
                Parcels = parcelpath;
            else
                Parcels = cifti_read(parcelpath);
                Parcels = Parcels.cdata;
            end
            ParcelCommunities =cifti_read(templatepath); % still use the Gordon colors for an unknown parcel?
            newCons.SortClus = zeros(length(Parcels),size(SortClus,2));
            for j = 1:size(SortClus,2)
                for i = 1:max(SortClus(:,j))
                     newCons.SortClus(any(Parcels==find(SortClus(:,j)==i)',2),j) = i;
                end
            end
            [CW,GenOrder] =assign_networks_by_template_cifti(newCons.SortClus,ParcelCommunities,0.1,'dice');%'dice'
        end
    case 4
         [CW,GenOrder] = makeCW(params,Nnets);
end

%% Re-Order Networks (vis,DMN,Mot,DAN,FPC,...)
% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=SortClus;foo(foo==0)=NaN;
SortClusRO=SortClus;
for j=1:length(GenOrder),SortClusRO(foo==GenOrder(j))=j;end
%% Check the colormaps are fine
figure;
imagesc;colormap(CWro.cMap);
colorbar;
title(sprintf('N = %i',size(CWro.cMap,1)))
pause(1);
close all;
