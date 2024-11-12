% this function was taken out from assign_network_colors for sorting the
% networks to a specific given order (e.g. anterior to posterior).
function SortClusRO = SortToGenOrder_HSB(SortClus,GenOrder)
foo=SortClus;foo(foo==0)=NaN;
SortClusRO=SortClus;
for j=1:length(GenOrder),SortClusRO(foo==GenOrder(j))=j;end
end