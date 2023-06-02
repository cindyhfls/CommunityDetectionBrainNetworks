function Gout=Fix_US_Grouping_HSB(Min,k)
%
% This function takes an ROI or voxel-wise InfoMap grouping matrix
% Min that is (Nroi,Nkden), selects a column to base the model off of, and
% then fills in grouping assignments for the ROIs in col k that are listed
% as 'UnSpecified' (defined as group==80).
% The algorithm steps to higher edge density (lower correlation) to test
% for grouping membership. If the group matches one that already exists,
% then the renaming of that ROI is complete. If not, it checks to see if
% the number of new elements in that group surpass the threshold (N=5 for
% ROI). If it does, then the new group is added.

%% Parameters
us=0; % Grouping label for UnSpecified group
th=5; % threshold for number of neighbors in group to keep

%% Choose initial starting column
Gout=squeeze(Min(:,k));

%% Find 'Unspecified' elements
toFix=find(Gout==us);
N2fix=length(toFix);

%% Redefine
if N2fix>0
New=[];
for j=1:N2fix
   posse=Min(toFix(j),:);           % Get entire row from matrix
   posse(posse==us)=[];             % Remove Junk
   if ~isempty(posse)
       newVal=mode(posse(:));
       if any(ismember(Gout,newVal))
           Gout(toFix(j))=newVal; % most frequent label in rest
       else
           New=cat(1,New,[j,newVal]);
       end
   end
end
if numel(New)<1, return, end

%% Test for new groups that are above size threshold
ux=unique(New(:,2));
if length(ux)==1
    counts=length(New(:,2));
else counts=hist(New(:,2),ux);
end

%% Add new groups if any above size threshold
if any(counts>=th)
    toRep=find(counts>=th);
    NtoRep=length(toRep);
for k=1:NtoRep
    idx=New(New(:,2)==ux(toRep(k)),1);
    Gout(toFix(idx))=ux(toRep(k));
end
end
ignore=ux(counts<th)';


%% Test remaining elements to see if they are ever in current group
toFix2=find(Gout==us);
N2fix=length(toFix2);
for j=1:N2fix
   posse=Min(toFix2(j),:); 
   posse=setdiff(posse,[us,ignore],'stable');
   if ~isempty(posse)
   Gout(toFix2(j))=posse(1);       
   end
end
end