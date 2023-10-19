function RGB2 = change_color_rgb_ncolor(RGB,ncolor)
    HSV = rgb2hsv(RGB);
    HSV = repelem(HSV,ncolor,1);
    HSV(:,2) = linspace(0.2,1,ncolor);
%     if HSV(1,3)>0.5
%         HSV(:,3) = linspace(0.5,1,ncolor);
%     end
    RGB2 = hsv2rgb(HSV);
end