function IM = make_template_from_parcel(Parcels,ParcelCommunities,ROIxyz)
% ROI assignments
% parcel_name = 'eLABE_Y2_prelim_05062023'
% [parcels_path] = Util.get_parcellation_path(parcel_name);
% Parcels = ft_read_cifti_mod(parcels_path);
% load(['/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/Parcels_',parcel_name,'.mat'],'ROIxyz');
% ParcelCommunities = ft_read_cifti_mod('/data/wheelock/data1/people/Cindy/BCP/Infomap/InfantTemplates/Kardan2022_communities.dlabel.nii');

%%
nParcels = length(setdiff(unique(Parcels.data),0));
net_id = NaN(nParcels,1);
for i = setdiff(unique(Parcels.data),0)'
    i
    net_id(i) = mode(ParcelCommunities.data(Parcels.data==i));
end

%%
clear IM
IM.name=['IM_temp'];
cMap = Plot.linspecer(length(unique(net_id)));
IM.cMap=cMap;
IM.Nets= string(1:size(cMap,1))';
IM.ROIxyz=ROIxyz;
Nroi = length(net_id);
IM.key=[[1:Nroi]',zeros(Nroi,1)];
[IM.key(:,2),IM.order]=sort(IM_Remove_Naming_Gaps_HSB(net_id));
IM.ROIxyz=IM.ROIxyz(IM.order,:);
IM=Org_IM_DVLR_HSB(IM);
end