function RGB2 = change_rgb_color(RGB,existing_colors)
if ~exist('existing_colors','var')||isempty(existing_colors)
    existing_colors = RGB;
end
% this function make a slightly different color than the input color
perceptual_thre = 20;
HSV = rgb2hsv(RGB);
tic
while true
    % Adjust Hue (considering wrap-around)
    HSV(1) = mod(HSV(1) + 0.01, 1);
    
    if HSV(2)<0.5 % saturation
        HSV(2) = min(HSV(2)+0.2,1);
    else
        HSV(2) = max(HSV(2)-0.2,0);
    end
    if HSV(3)<0.5 % brightness
        HSV(3) = min(HSV(3)+0.2,1);
    else
        HSV(3) = max(HSV(3)-0.2,0);
    end
    
    RGB2 = hsv2rgb(HSV);
    
    % convert to CIELAB
    LAB2 = rgb2lab(RGB2);
    
    % compute perceptual difference
    perceptual_diffs = zeros(size(existing_colors,1),1);
    for i = 1:size(existing_colors,1)
        LAB1 = rgb2lab(existing_colors(i,:));
        perceptual_diffs(i) = sqrt(sum((LAB1-LAB2).^2));
    end
    t = toc;
    if all(perceptual_diffs>perceptual_thre)||(t>200)
        break
    end
end

end