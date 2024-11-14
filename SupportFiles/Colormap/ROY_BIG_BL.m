function map = ROY_BIG_BL(n)

% custom function to reproduce the ROY_BIG_BL in connectome workbench
if nargin < 1
    n = size(get(gcf, 'Colormap'), 1);
end

values = load('cmap_ROY_BIG_BL.mat');
values = values.cmap_ROY_BIG_BL;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');

end
