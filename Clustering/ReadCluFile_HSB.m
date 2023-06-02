function M=ReadCluFile_HSB(fn)
%
% This function reads in a *.clu file output from infomap (current version
% 0.18.9 as of 9/20/2016).
% The output M is organized as [node idx, module label, total flow].
% Flow is a metric used by the map equation to determine the module
% structure. See http://www.mapequation.org/code.html for more info.
% The output M is sorted to be in the same ordering as the input nodes in
% the pajek file.


fid=fopen(fn,'r');
tsc=textscan(fid,'%f %f %f','HeaderLines',2);
fclose(fid);
M=[double(tsc{1}),double(tsc{2}),double(tsc{3})];
M=sortrows(M,1);