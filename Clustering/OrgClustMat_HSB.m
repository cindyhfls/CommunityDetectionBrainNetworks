function matOUTa=OrgClustMat_HSB(matIN,killTH)
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


%% Set grouping value == 0 if less than threshold membership
matIN = remove_singleton(matIN,killTH);

%% Hungarian Matching
for j = 2:size(matIN ,2)
    matIN (:,j) = pair_labeling(matIN (:,j-1),matIN (:,j));
end
matOUT = matIN;
%% Initialize (removed)
% maxclr=max(matIN(:));
% matKey=matIN;
% matKey(matKey==junk)=NaN;
% matOUT=zeros(Nroi,Nkden,'single');
% matOUT(:,1)=matIN(:,1); % 1st col is key
%% Cycle through other columns (removed)
% for j=2:Nkden
%     key=squeeze(matOUT(:,(j-1)));   % earlier column values
%     KeyVals=setdiff(unique(key),0); % update 2023.09.28 JCT
%     Nold=length(KeyVals);
%     new=squeeze(matKey(:,j));        % to-re-organize column values
%     NewVals=unique(new(~isnan(new)));% update 2023.09.28 JCT
%     Nnew=length(NewVals);
%     
%     % Make a matrix of # overlaps
%     ol=zeros(Nold,Nnew,'single');
%     for l=1:Nold
%         for k=1:Nnew
%             ol(l,k)=sum((key==KeyVals(l)).*(new==NewVals(k)));
%         end
%     end
%     
%     
%     [olp(:,1),olp(:,2),olp(:,3)]=find(ol); % Get sorted overlaps
%     cc=sortrows(olp,-3);
% 
%     u1=zeros(Nold,1);       % Assign r2 to r1 values when r1 unused
%     u2=zeros(Nnew,1);
%     while ~isempty(cc)
%     
%      % if both old and new colors in question haven't been assigned
%         if (u1(cc(1,1))==0 && u2(cc(1,2))==0)
%             
%             u1(cc(1,1))=1;  % now they're assigned
%             u2(cc(1,2))=1;
%             
%             % old clrs with this value will be reassigned
%             exitnum=NewVals(cc(1,2));
%             % to this value from the newclrs
%             insertnum=KeyVals(cc(1,1));
%             % and here they are assigned
%             matOUT((matKey(:,j)==exitnum),j)=insertnum;
%             
%             % this deletes all entries of this oldclr from the unique list
%             cc=cc(~(cc(:,2)==cc(1,2)),:);
%     
%     % if an oldclr can't find parthers in newclrs (the newclr has
%     % already been given away, probably)
%         elseif (u1(cc(1,1))==1 && u2(cc(1,2))==0)
%             
%             % now this color is used from oldclrs
%             u2(cc(1,2))=1;
%             % and we will be replacing oldclrs with this value
%             exitnum=NewVals(cc(1,2));
%             % with the highest number possible
%             matOUT((matKey(:,j)==exitnum),j)=maxclr;
%             % and then we'll bump the highest number possible up one more
%             maxclr=maxclr+1;
%             % and now we'll delete entries of this oldclr from the list
%             cc=cc(~(cc(:,2)==cc(1,2)),:);
%         
%         % the case before this case should take care of all scenarios, but
%         % this is a catch in case somehow a NewVals didn't find a partner in
%         % KeyVals, and the oldclr has already been used
%         elseif (u1(cc(1,1))==0 && u2(cc(1,2))==1)
%             disp newold! how did i get in here?
%             cc=cc(~(cc(:,2)==cc(1,2)),:);
%         end
%     clear exitnum insertnum;
%     end
%     clear olp
% end
% matOUT(~isfinite(matOUT))=0;
%%
% now drop all values to the minimum possible
vals=unique(matOUT);
matOUTa=zeros(Nroi,Nkden);
for j=1:size(vals,1)
    matOUTa(matOUT==vals(j))=j;
end
if any(vals==0)
    matOUTa=matOUTa-1;% keep junk as zero
end
matOUTa=single(matOUTa);