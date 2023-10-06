% Example
% Parcels =ft_read_cifti_mod('/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/eLABE_Y2_N113_atleast600frames_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg_global_edgethresh_0.75_nogap_minsize_15_relabelled.dlabel.nii');
% Parcels = Parcels.data;
% NewParcels = ft_read_cifti_mod('/data/wheelock/data1/parcellations/333parcels/Parcels_LR.dtseries.nii');
% NewParcels = NewParcels.data;
% load('/data/wheelock/data1/people/Cindy/BCP/Infomap/parcel-wise/eLABE_Y2_N113/eLABE_Y2_prelim_072023_0.75/230927/IM_Infomap_eLABE_Y2_N113_low0.001_step0.001_high0.100_xdist20.mat_Consesus_model_2.mat')


function [IMnew] = convert_key_parcel_to_parcel(IM, Parcels, NewParcels)
[~,sortorder] = sort(IM.order);
key = IM.key(sortorder,2);
assn = NaN(size(Parcels));
for k = unique(key)'
    assn(any(Parcels==find((key==k))',2))=k;
end
newkey = NaN(size(setdiff(unique(NewParcels),0)'))';
for i = setdiff(unique(NewParcels),0)'
    tmp = assn(NewParcels==i);
    newkey(i) = mode(setdiff(tmp,[NaN,0]));
end

newkey(isnan(newkey)) = Inf;
IMnew.name=['IM_temp'];
keepnets = setdiff(unique(newkey),Inf);
keepnets = keepnets(~isnan(keepnets));
IMnew.cMap=IM.cMap(keepnets,:);
IMnew.Nets = IM.Nets(keepnets);
Nroi = length(newkey);
IMnew.key=[[1:Nroi]',zeros(Nroi,1)];
[IMnew.key(:,2),IMnew.order]=sort(IM_Remove_Naming_Gaps_HSB(newkey));
if sum(newkey==Inf)>0
    IMnew.key(IMnew.key(:,2)==max(IMnew.key(:,2)),2)=0;
end
end
