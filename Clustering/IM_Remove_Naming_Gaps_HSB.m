function [foo,pairing]=IM_Remove_Naming_Gaps_HSB(foob)

% This function removes gaps in group numbering lables. The output is
% numbered in order from 1 to the number of networks.
foo=foob;
Nets=unique(foob(:));
Nnets=length(Nets);
for j=1:Nnets
    foo(foob==Nets(j))=j;
end
pairing=[Nets',1:Nnets];