function cmap = make_sequential_colors(start_color,end_color,n)
% Number of colors in the colormap
if nargin<3
    n = 256;
end

% Create red, green, and blue components
r = linspace(start_color(1), end_color(1), n);
g = linspace(start_color(2), end_color(2), n);
b = linspace(start_color(3), end_color(3), n);

% Combine them into a colormap
cmap = [r' g' b'];

end