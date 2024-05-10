%% Viewing
iNet = 2
Edit_NetworkColors(stats.SortClusRO,CWro,iNet,Parcels);
%% Editing
customcolor = [1,0,0];% fill in if you want to change it
% customcolor = distinguishable_colors(1,CWro.cMap); % set color to one
CWro2=Edit_NetworkColors(stats.SortClus,CWro,iNet,Parcels,customcolor);