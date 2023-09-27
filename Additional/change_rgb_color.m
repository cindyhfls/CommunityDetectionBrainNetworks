function RGB2 = change_rgb_color(RGB)
% make a slightly different color

HSV = rgb2hsv(RGB);
% if HSV(1)<0.5 % hue
%     HSV(1) = HSV(1)+0.1;
% else
%     HSV(1) = HSV(1)-0.1;
% end
if HSV(2)<0.5 % saturation
    HSV(2) = HSV(2)+0.2;
else
    HSV(2) = HSV(2)-0.2;
end
if HSV(3)<0.5 % brightness
    HSV(3) = HSV(3)+0.2; 
else
    HSV(3) = HSV(3)-0.2;
end

HSV(HSV > 1) = 1;  % Limit values
RGB2 = hsv2rgb(HSV);

end