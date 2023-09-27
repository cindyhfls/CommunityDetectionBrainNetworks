function [foo, CWro,Cons] = assign_network_colors(Cons,nameoption,templatepath)
%%
if ~exist('templatepath','var')||isempty(templatepath)
    templatepath = 'IM_Gordon_2014_333_Parcels.mat';
end

%% Label and Color Brain Networks identified by Infomap

% Find the unique networks identified by Infomap
Nets=unique(Cons.SortCons(:));
Nnets=length(Nets);

% nameoption = 3
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
        %%
%         parcel_name = 'Gordon'%'eLABE_Y2_prelim_05062023'
%         [parcels_path] = Util.get_parcellation_path(parcel_name);
%         Parcels = ft_read_cifti_mod(parcels_path);
%         load(['Parcels_',parcel_name,'.mat'],'ROIxyz');
%         ParcelCommunities = ft_read_cifti_mod('Kardan2022_communities.dlabel.nii');
%         ParcelCommunities = ft_read_cifti_mod('Parcel_Communities_cortexonly_12Networks.dlabel.nii')
%         colortemplate.IM = make_template_from_parcel(Parcels,ParcelCommunities,ROIxyz);
        colortemplate =load(templatepath);
%         template = load('/data/wheelock/data1/parcellations/IM/Kardan_2022_DCN/IM_11_BCP94.mat');
        % to-do: convert template from the cifti(see how I made the
        % Wange Network parcels)
        [CW,GenOrder,MIn] = assign_Infomap_networks_by_template(Cons,colortemplate,0.1,'dice');
end

%% Re-Order Networks (vis,DMN,Mot,DAN,FPC,...)
% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=Cons.SortCons;foo(foo==0)=NaN;
Cons.SortConsRO=Cons.SortCons;
for j=1:length(GenOrder),Cons.SortConsRO(foo==GenOrder(j))=j;end
foo=Cons.SortConsRO;
