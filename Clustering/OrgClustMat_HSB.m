function matOUTa=OrgClustMat_HSB(matIN,killTH,reverse)
%
% This function organizes a matrix representing clustering/sorting
% assignments. The input matrix rows are ROIs sorted within column. The
% numbers representing groups within each column are not related across
% column. This matrix re-organizes by 1st column so that there is
% consistency in labels (numbers). It is assumed that matIN is 2D.
% The first column is the master key. Thus, this code could be used to
% re-sort a matrix based on a 'key' sorting input as its 1st column.
% In the output matrix, ZERO = junk, aka UnSpecified


%% Set parameters
junk=0;
[Nroi,Nkden]=size(matIN);
if ~exist('killTH','var')
    if Nroi<1000 % ROI case    
        killTH=5;
    else
        killTH=100; % MVW/node cases
    end
end
if ~exist('reverse','var')
    reverse =0;
end
matOUT = remove_singleton(matIN,killTH);
[matOUT] =  regularize_HSB(matOUT,reverse);

% [matOUT] =  regularize_HSB(matIN,reverse);
% %% Hungarian Matching (somehow this sometimes give completely
% non-overlapping communities the same index)
% if reverse
%     for j = size(matIN,2)-1:-1:1
%         matIN (:,j) = pair_labeling(matIN (:,j+1),matIN (:,j));
%     end
% else
%     for j = 2:size(matIN ,2)
%         matIN (:,j) = pair_labeling(matIN (:,j-1),matIN (:,j));
%     end
% end
% matOUT = matIN;
%%

% % 1) merge the ones with high overlap but are non-continous
% G1=setdiff(unique(matOUT(:)),0);
% for j=1:length(G1)
%     [G1(j,2),G1(j,3)]=max(sum(matOUT==G1(j,1)));
% end
% repnets = cell2mat(arrayfun(@(ii)matOUT(:,G1(ii,3))==G1(ii,1),[1:length(G1)]','UniformOutput',false)');
% 
% % calculate overlap
% ovlp = zeros(size(repnets,2));
% for i = 1:size(repnets,2)
%     for j = i+1:size(repnets,2)
%         ovlp(i,j) = dice(repnets(:,i),repnets(:,j));
%     end
% end
% [idxi,idxj]= find(ovlp>0.8); % arbitrary threshold to say that dice>0.5 are too similar to be considered as two networks
% 
% for ii = 1:length(idxi)
%     lowid = min([idxi(ii),idxj(ii)]);
%     highid = max([idxi(ii),idxj(ii)]);
%     if isempty(intersect(find(sum(matOUT==G1(lowid,1))~=0),find(sum(matOUT==G1(highid,1))~=0)))% if they don't appear in the same solutions ever
%         fprintf('Merge network %i and %i\n',G1(lowid,1),G1(highid,1));
%         matOUT(matOUT==G1(highid,1)) = G1(lowid,1); % reassign the higher number to the lower number
%     end
% end

% 2) remove networks that only appear in one threshold
% appearance = arrayfun(@(ii)max(sum(matOUT==ii,2)),setdiff(unique(matOUT(:)),0));
% for k = find(appearance==1)'
%     matOUT(matOUT==k) = 0;
% end

% 3) remove networks with lower than killTH number of nodes
% matOUT = remove_singleton(matOUT,killTH);

%% now drop all values to the minimum possible
vals=unique(matOUT);
matOUTa=zeros(Nroi,Nkden);
for j=1:size(vals,1)
    matOUTa(matOUT==vals(j))=j;
end
if any(vals==0)
    matOUTa=matOUTa-1;% keep junk as zero
end
matOUTa=single(matOUTa);