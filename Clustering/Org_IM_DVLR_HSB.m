function IM=Org_IM_DVLR_HSB(IM)
%
% This function organizes the ROIs in an IM structure within networks in
% order of dorsal to ventral within left then right.
%

%% do it
Nnets=size(IM.Nets,1);
for j=1:Nnets
    % General network ROIs
    idxNet=find(IM.key(:,2)==j);
    % L ROIs
    idxL=find(IM.ROIxyz(idxNet,1)<0);
    [~,idxLdvI]=sort(IM.ROIxyz(idxNet(idxL),3));
    % R ROIs
    idxR=find(IM.ROIxyz(idxNet,1)>=0);
    [~,idxRdvI]=sort(IM.ROIxyz(idxNet(idxR),3));
    
    idxNetRO=[idxL(idxLdvI);idxR(idxRdvI)];
    
    IM.ROIxyz(idxNet,:)=IM.ROIxyz(idxNet(idxNetRO),:);
    IM.order(idxNet)=IM.order(idxNet(idxNetRO));
end