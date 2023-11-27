function [clusters,codelengthTwoLevels]=ReadCluFile_HSB(fn,version)
%
% This function reads in a *.clu file output from infomap (current version
% 0.18.9 as of 9/20/2016).
% The output M is organized as [node idx, module label, total flow].
% Flow is a metric used by the map equation to determine the module
% structure. See http://www.mapequation.org/code.html for more info.
% The output M is sorted to be in the same ordering as the input nodes in
% the pajek file.
% 2023.11.10 adding support for the original 0.15.7 version

if ~exist('version','var')||isempty(version)
    version = '0.18.9'; % default version
end

clusters = NaN; codelengthTwoLevels = NaN;

if strcmp(version,'0.18.9')
    fid=fopen(fn,'r');
    l = fgetl(fid);
    tsc=textscan(fid,'%f %f %f','HeaderLines',1);
    fclose(fid);
    M=[double(tsc{1}),double(tsc{2}),double(tsc{3})];
    M=sortrows(M,1);
    clusters=M(:,2);
    % Extract codelength in two levels
    codelengthTwoLevels_match = regexp(l, 'codelength (\d+\.\d+) in 2 levels', 'tokens');
    codelengthTwoLevels = str2double(codelengthTwoLevels_match{1});
elseif strcmp(version,'0.15.7') % J.Powers version
    fid = fopen(fn,'r');
    tsc= textscan(fid,'%d','HeaderLines',1);
    clusters = tsc{1};
else
    error('currently only supporting reading from version 0.18.9 or 0.15.7, please edit ReadCluFile_HSB.m to incorporate reading from the .clu file from Infomap')
end