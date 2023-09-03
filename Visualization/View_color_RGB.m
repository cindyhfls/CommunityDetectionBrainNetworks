% 2022.08.23
% adapted from https://www.mathworks.com/matlabcentral/answers/409780-hi-guys-how-can-i-have-display-a-color-from-rgb-coordinates
function View_color_RGB(clr,varargin)

name = '';
if ~isempty(varargin)
    name = varargin{1};
end

if max(clr)>1
    clr = clr/255;
end
x = [0 1 1 0] ; y = [0 0 1 1] ;

figure;
fill(x,y,clr)
title(name)
end