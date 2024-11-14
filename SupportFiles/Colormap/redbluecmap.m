function map = f_redbluecmap(n)
if nargin < 1
    n = size(get(gcf, 'Colormap'), 1);
end

c = [103    0       31;
    178     24      43;
    214     96      77;
    244     165     130;
    253     219     199;
    247     247     247;
    209     229     240;
    146     197     222;
    67      147     195;
    33      102     172;
    5       48      97];
values = flipud(c/255);

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');

end
