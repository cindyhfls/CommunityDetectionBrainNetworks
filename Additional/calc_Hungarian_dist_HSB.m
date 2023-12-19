function Hungarian_distance = calc_Hungarian_dist_HSB(ref_labels,input_labels)
% Adapted from CBIG_HungarianClusterMatch.m Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

assert(isvector(ref_labels) && isvector(input_labels))

n= length(ref_labels);

num_ref_labels = max(ref_labels);
num_input_labels = max(input_labels);

% Build matching matrix
mat = zeros(num_input_labels, num_ref_labels);
for i = 1:num_ref_labels
    for j = 1:num_input_labels
        mat(j, i) = -sum(double(ref_labels(:) == i) .* double(input_labels(:) == j));
    end
end

[assign, cost] = munkres(mat);
Hungarian_distance = 1-(-cost/n);

end