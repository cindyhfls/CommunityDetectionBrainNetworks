function[] = mat_plot(A)

set(gcf, 'color', 'w');
colormap(jet(256));
imagesc(A);
colorbar;
axis image;