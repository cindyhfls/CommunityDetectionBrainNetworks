function l = hline(y0, varargin)
%In the newer MATLAB versions you can use yline(y0).
if nargin < 1
  y0 = 0;
end
a = axis;
l = line([a(1) a(2)], [y0 y0]);
set(l, 'LineStyle', '--', 'LineWidth', 0.5, 'color', 'k');

if length(varargin) > 0
  set(l, varargin{:});
end
end
